'''
Classes used to parse and store VCF data
'''
import pysam, copy
import pandas as pd
from gumpy import GeneticVariation


class VCFRecord(object):
    '''
    Class for VCF records
    Instance variables:
        `chrom` (str): Name of the sample.
        `pos` (int): Genome index for the change.
        `ref` (str): Reference value for the nucleotide.
        `alts` (tuple(str)): Alternative calls. Tuple can contain single values and indels.
        `qual` (None): Values for quality (values other than None have yet to be found in testing files...).
        `filter` (str): Whether this record has pass the filter.
        `info` (dict): Dictionary of key->value for the info fields.
        `values` (dict): Dictionary of key->value for the values. Usually this is the FORMAT field names with their corresponding values.
        `is_filter_pass` (bool): does the filter column contain PASS?
        `call1` (int): the index of the first call
        `call2` (int): the index of the second call
        `is_reference` (bool): is the call for the reference?
        `is_null` (bool): is the call a null call?
        `is_heterozygous` (bool): is it a is_heterozygous call i.e. call1!=call2?
        `is_alt` (bool): or, is the call for a single specified alt
    '''
    def __init__(self, *args, **kwargs):
        '''Constructor, pulls the data out of the pysam object
        into a more useable format

        Args:
            record (pysam.libcbcf.VariantRecord): The record object
            sample_num (str) : Name of the sample to consider. Used for possible cases
                                where there is more than 1 sample per record
        '''
        if len(args) != 2:
            #Rebuilding...
            assert "reloading" in kwargs.keys(), "Incorrect arguments given."
            allowed_kwargs = ['chrom', 'pos', 'ref', 'alts', 'qual', 'filter', 'info', 'values']
            seen = set()
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
                    seen.add(key)
            for key in set(allowed_kwargs).difference(seen):
                #Default values to None if the kwarg has not been passed
                setattr(self, key, None)
            return
        else:
            record = args[0]
            sample = args[1]

        assert len(record.samples.keys())==1, 'only supporting single samples per row at present!'

        assert 'GT' in record.samples[sample].keys(), 'require GT in FORMAT column to parse genotype'

        #Save some of the easier to get to attributes
        self.chrom = record.chrom
        self.contig = record.contig
        self.pos = record.pos
        self.ref = record.ref.lower()
        if isinstance(record.alts,tuple):
            self.alts = tuple([i.lower() for i in record.alts])
        else:
            self.alts=record.alts
        self.qual = record.qual

        #Get the filter attribute value
        assert len(record.filter.items()) >= 0, "A record has more than 1 filter set!"
        if len(record.filter.items()) == 0:
            self.filter = None
            self.is_filter_pass=False
        else:
            self.filter = record.filter.items()[0][0]
            self.is_filter_pass=True if record.filter.items()[0][0] == 'PASS' else False

        #Get the info field
        self.info = {}
        for (key, value) in record.info.items():
            self.info[key] = value

        #Get the values
        self.values = {}
        for (key, item) in record.samples[sample].items():
            # incorporate the logic from the old Genotype class here
            if key == 'GT':
                call1,call2=item
                self.call1 = call1 if call1 is not None else -1
                self.call2 = self.call1 if call2 is None else call2
                self.is_reference = True if call1==0 and (call1 == call2 or call2 == -1) else False
                self.is_heterozygous = True if call1!=call2 else False
                self.is_null = True if set([self.call1, self.call2]) == {-1} else False
                self.is_alt = True if not self.is_reference and not self.is_heterozygous and not self.is_null else False

            #Due to how pysam reads floats, there are some erroneously long dps, so round
            elif isinstance(item,float):
                item = round(item, 3)
            self.values[key] = item

    def __repr__(self):
        '''Pretty print the record
        '''
        s = self.chrom + "\t"
        s += str(self.pos) + "\t"
        s += self.ref + "\t"
        s += str(self.alts) + "\t"
        if self.qual == 'None':
            s+='.\t'
        else:
            s+=str(self.qual)+'\t'
        if self.filter == None:
            s+='.\t'
        else:
            s+=self.filter+'\t'
        for val in self.values.keys():
            s += str(val)+":"
        s = s[:-1] + "\t"
        for val in self.values.values():
            s += str(val)+":"
        s = s[:-1]
        s += "\n"
        return s

class VCFFile(object):
    '''
    Class to instanciate a variant file (VCF)
    Used to apply a VCF file to a genome
    Instance variables:
        `vcf_version` (tuple(int)): Tuple of ints to show the VCF version of the file. e.g 3.2 would be (3, 2).
        `contig_lengths` (dict): Dictionary mapping contig_name->length for all defined contigs.
        `formats` (dict): Dictionary mapping format_name->dict(description, id, type).
        `records` (list(VCFRecord)): List of VCFRecord objects for each record within the file.
        `changes` (dict): Dictionary mapping array_index->(call, [(n_reads, alt_call)])
    '''
    def __init__(self, *args, **kwargs):
        '''
        Constructor. Reads the VCF file and sets up values
        Args:
            filename (str) : The name of the VCF file
        '''
        if len(args) != 1:
            #Rebuilding...
            assert "reloading" in kwargs.keys(), "Incorrect arguments given. Only give a filename."
            allowed_kwargs = ['vcf_version', 'contig_lengths', 'formats', 'records', 'changes', 'ignore_filter', 'format_fields_min_thresholds', 'bypass_reference_calls']
            seen = set()
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
                    seen.add(key)
            for key in set(allowed_kwargs).difference(seen):
                #Default values to None if the kwarg has not been passed
                setattr(self, key, None)
            return
        else:
            filename = args[0]

        self.ignore_filter=kwargs.get('ignore_filter',False)
        self.bypass_reference_calls=kwargs.get('bypass_reference_calls',False)
        self.format_fields_min_thresholds=kwargs.get('format_fields_min_thresholds',None)
        self.filename=filename

        #Use pysam to parse the VCF
        vcf = pysam.VariantFile(self.filename)

        #Get some basic metadata
        self.vcf_version = vcf.version

        #Get the contig lengths from the header
        self.contig_lengths = {}
        for name in list(vcf.header.contigs):
            self.contig_lengths[name] = vcf.header.contigs[name].length

        #Get the formats
        self.format_fields_metadata = {}
        for format in vcf.header.formats.keys():
            description = vcf.header.formats[format].description
            id = vcf.header.formats[format].id
            f_type = vcf.header.formats[format].type
            self.format_fields_metadata[format] = {
                "description" : description,
                "id": id,
                "type": f_type
            }

        if isinstance(self.format_fields_min_thresholds,dict):
            assert set(self.format_fields_min_thresholds.keys()).issubset(set(self.format_fields_metadata.keys())), "field to threshold on not found in the FORMAT column of the vcf!"

        #Get the records
        self.records = []
        for record in list(vcf):
            for sample in record.samples.keys():
                self.records.append(VCFRecord(record, sample))

        #Ensure that only a single record exists for each position specified
        assert len(self.records) == len(set([record.pos for record in self.records])), "There must only be 1 record per position!"

        self.__find_calls()

    def __repr__(self):
        '''Pretty print the VCF file

        Returns:
            str: String summarising the VCF file
        '''
        output='VCF variant file, version '+''.join(str(i)+"." for i in self.vcf_version)[:-1]+'\n'
        output+=self.filename+'\n'
        output+=str(len(self.records))+' records'+'\n'
        output+='FORMAT columns: '+', '.join(i for i in sorted(list(self.format_fields_metadata.keys())))+'\n'+'\n'
        if len(self.records)>3:
            output += str(self.records[0])
            output += str(self.records[1])
            output += str(self.records[2])
            output += "...\n"
            output += str(self.records[-1])
        else:
            for record in self.records:
                output += str(record)+'\n'
        return output

    def __find_calls(self):
        '''Function to find changes within the genome based on the variant file
        '''
        self.calls = {}
        for record in self.records:

            # VCF files are 1 indexed but keep for now
            index = copy.deepcopy(record.pos)

            # if we've asked, bypass (for speed) if this is a ref call
            if self.bypass_reference_calls and record.is_reference:
                continue

            # bypass filter fails , unless we have asked to ignore filter calls
            if not self.ignore_filter and not record.is_filter_pass:
                continue

            # only proceed if a dictionary has been passed (otherwise defaults to None)
            proceed=True
            if isinstance(self.format_fields_min_thresholds,dict):
                # ok to just do since we've already check in the constructor that these fields exist in the VCF
                for i in self.format_fields_min_thresholds:
                    proceed = proceed and record.values[i] >= self.format_fields_min_thresholds[i]
            if not proceed:
                continue

            if record.is_heterozygous:
                variant='z'*len(record.ref)
                variant_type='het'
            elif record.is_alt:
                variant=record.alts[record.call1-1]
                variant_type='snp'
            elif record.is_null:
                variant='x'*len(record.ref)
                variant_type='null'
            elif record.is_reference:
                variant=record.ref
                variant_type='ref'

            if len(record.ref)==len(variant):
                for counter,(before,after) in enumerate(zip(record.ref,variant)):

                    # only make a change if the ALT is different to the REF
                    if before!=after:
                        metadata={}
                        metadata['type']=variant_type
                        metadata['call']=after
                        metadata['ref']=before
                        metadata['pos']=counter
                        vcf_info={}
                        vcf_info=copy.deepcopy(record.values)
                        vcf_info['REF']=record.ref
                        vcf_info['ALTS']=record.alts
                        metadata['original_vcf_row']=vcf_info
                        self.calls[index+counter]=metadata

            else:
                mutations = self._simplify_call(record.ref, variant)
                for (p, type_, bases) in mutations:
                    # p = max(p - 1, 0)
                    if type_ in ["ins", "del"]:
                        indel_length = len(bases)
                        # if type_ == "del":
                        #     indel_length *= -1
                        metadata = {}
                        metadata['type'] = 'indel'
                        metadata['call'] = (type_,bases)
                        # metadata['call'] = (type_, indel_length)
                        metadata['ref']=record.ref[0]
                        metadata['pos']=p
                        vcf_info={}
                        vcf_info=copy.deepcopy(record.values)
                        vcf_info['REF']=record.ref
                        vcf_info['ALTS']=record.alts
                        metadata['original_vcf_row']=vcf_info
                        self.calls[index+p]=metadata
                    else:
                        metadata = {}
                        metadata['type'] = variant_type
                        metadata['call'] = bases[1]
                        metadata['ref'] = bases[0]
                        metadata['pos'] = p
                        vcf_info = {}
                        vcf_info = copy.deepcopy(record.values)
                        vcf_info['REF'] = record.ref
                        vcf_info['ALTS'] = record.alts
                        metadata['original_vcf_row'] = vcf_info
                        self.calls[index+p] = metadata

    def _simplify_call(self, ref, alt):
        '''Find where in the sequence the indel was, and the values.
        Based on finding the indel position at which there is the least SNPs

        Args:
            ref (str): Reference bases. Should match reference bases at this point
            alt (str): Alt bases.
        Returns:
            [(int, str, str)]: Returns a list of tuples of (pos, one of ['ins','del','snp'], indel_bases or (ref, alt) for SNPs)
        '''
        def snp_number(ref, alt):
            '''Count the number of SNPs between 2 sequences

            Args:
                ref (str): Reference bases
                alt (str): Alt bases

            Returns:
                int: Number of SNPs between ref and alt
            '''
            snps = 0
            for (a, b) in zip(ref, alt):
                if a is not None and b is not None and a != b:
                    snps += 1
            return snps

        '''The process for finding the positions for indels are almost identical
        as the process for finding a del can be interpreted as finding an ins with ref and alt reversed.
        The approach here is to use a sliding window to find the position of the indel where the number of SNPs is lowest.
        Assumes that there is only a single indel between the ref and the alt - there may be cases where this does not work
            these will just produce large amounts of SNPs... Could be adapted to check for multi-indel but this will scale
            exponentially the number of versions checked.
        '''
        if len(ref) > len(alt):
            #Del
            length = len(ref) - len(alt)
            x = ref
            y = alt
            indel = "del"
        elif len(ref) < len(alt):
            #Ins
            length = len(alt) - len(ref)
            x = alt
            y = ref
            indel = "ins"
        start = 0
        current = None
        current_snps = 999
        mutations = []
        for i in range(len(y)+1):
            y1 = [y[a] for a in range(i)] + [None for a in range(length)] + [y[a] for a in range(i, len(y))]
            if snp_number(y1, x) <= current_snps:
                current = y1
                current_snps = snp_number(y1, x)
                start = i
        seq = [x[i] for i in range(len(current)) if current[i] is None]
        #Add the indel
        # if indel == "ins":
            #Ins after index, del at index so adjust ins
            # start -= 1
        mutations.append((start, indel, ''.join(seq)))
        #Check for SNPs and add those
        for (i, (a, b)) in enumerate(zip(x, current)):
            if a is not None and b is not None and a != b:
                if indel == "ins":
                    # PWF was (i-1) but corrected, prob
                    mutations.append((i, "snp", (b, a)))
                else:
                    mutations.append((i, "snp", (a, b)))
        return mutations


    def to_df(self):
        '''Convert the VCFRecord to a pandas DataFrame.
        Metadata is stored in the `attrs` attribute of the DataFrame which may break with some operations
        (but pandas does not currently have a robust method for metadata storage...)

        Returns:
            pandas.DataFrame: DataFrame containing all of the information from the VCF file
        '''
        meta_data = {
            "vcf_version": self.vcf_version,
            "contig_lengths": self.contig_lengths,
            "formats": self.format_fields_metadata
        }

        chroms = []
        pos = []
        refs = []
        alts = []
        qual = []
        infos = []
        values = {}
        for record in self.records:
            chroms.append(record.chrom)
            pos.append(record.pos)
            refs.append(record.ref)
            alts.append(record.alts)
            qual.append(record.qual)
            infos.append(record.info)
            for key in record.values:
                if values.get(key) is None:
                    values[key] = [record.values[key]]
                else:
                    values[key].append(record.values[key])
        df = pd.DataFrame({
            "CHROM": chroms, "POS": pos,
            "REF": refs, "ALTS": alts,
            "QUAL": qual, "INFO": infos,
            **values
            })
        df.attrs = meta_data
        return df

    def interpret(self, genome):
        '''Takes in a Genome object, returning an object detailing the full differences caused by this VCF including amino acid differences

        Args:
            genome (gumpy.Genome): A reference Genome object
        '''
        return GeneticVariation(self, genome)
