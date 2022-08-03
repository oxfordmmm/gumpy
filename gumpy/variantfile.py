'''
Classes used to parse and store VCF data
'''
import copy
import pathlib
from collections import defaultdict

import numpy
import pandas
import pysam


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
        '''Constructor for the VCFRecord object.

        Parses the supplied pysam object and presents in a more Pythonic format

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
        `filename` (str): path to the VCF file
        `vcf_version` (tuple(int)): Tuple of ints to show the VCF version of the file. e.g 3.2 would be (3, 2).
        `contig_lengths` (dict): Dictionary mapping contig_name->length for all defined contigs.
        `format_fields_description` (dict): Dictionary mapping format_name->dict(description, id, type).
        `records` (list(VCFRecord)): List of VCFRecord objects for each record within the file.
        `calls` (dict): Dict of definite calls made in the VCF file, after any additional filtering has been applied
        `ignore_filter` (bool): whether to ignore the FILTER in the VCF file
        `format_fields_min_thresholds` (dict): dictionary specifying minimum thresholds to be applied to fields in the FORMAT field e.g. {'GTCONF':5}
        `variants` (numpy.array): Numpy array of the detected variants in the VCF file
        `nucleotide_index` (numpy.array): Array of genome indices which are affected by the VCF
        `ref_nucleotides` (numpy.array): Array of REF bases
        `alt_nucleotides` (numpy.array): Array of ALT bases
        `indel_length` (numpy.array): Array of lengths of insertions (+ve) or deletions (-ve) at each site
        `is_snp` (numpy.array): Array to act as a mask for `nucleotide_index` to show which are SNPs
        `is_het` (numpy.array): Array to act as a mask for `nucleotide_index` to show which are heterozygous calls
        `is_null` (numpy.array): Array to act as a mask for `nucleotide_index` to show which are null calls
        `is_indel` (numpy.array): Array to act as a mask for `nucleotide_index` to show which are indel calls
        `snp_distance` (int): SNP distance caused by the VCF
    '''
    def __init__(self, *args, **kwargs):
        '''
        Constructor for the VCFFile object.

        Parses the VCF file using pysam.

        Args:
            filename (str) : The name of the VCF file
            ignore_filter (bool): If True, ignore the FILTER column in the VCF file. Default is False.
            bypass_reference_calls (bool): If True, skip any row in the VCF (and therefore do not record)  which calls reference (i.e. 0/0). Default is False.
            format_fields_min_thresholds (dict, optional): Dict of field name in the FORMAT column and a minimum threshold to apply e.g. {'DP':5}
        '''
        if len(args) != 1:
            #Rebuilding...
            assert "reloading" in kwargs.keys(), "Incorrect arguments given. Only give a filename."
            allowed_kwargs = ['ignore_filter', 'format_fields_min_thresholds', 'bypass_reference_calls']
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
        assert isinstance(self.ignore_filter,bool)

        self.bypass_reference_calls=kwargs.get('bypass_reference_calls',False)
        assert isinstance(self.bypass_reference_calls,bool)

        self.format_fields_min_thresholds=kwargs.get('format_fields_min_thresholds',None)
        if self.format_fields_min_thresholds is not None:
            assert isinstance(self.format_fields_min_thresholds,dict)

        self.filename=filename
        assert isinstance(self.filename,str)
        assert pathlib.Path(self.filename).is_file()

        #Use pysam to parse the VCF
        #Pylint doesn't think there is a VariantFile method but it works...
        vcf = pysam.VariantFile(self.filename)

        #Get some basic metadata
        self.vcf_version = vcf.version

        #Get the contig lengths from the header
        self.contig_lengths = {}
        for name in list(vcf.header.contigs):
            self.contig_lengths[name] = vcf.header.contigs[name].length

        #Get the formats
        self.format_fields_metadata = {}
        for format_ in vcf.header.formats.keys():
            description = vcf.header.formats[format_].description
            id_ = vcf.header.formats[format_].id
            f_type = vcf.header.formats[format_].type
            self.format_fields_metadata[format_] = {
                "description" : description,
                "id": id_,
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
        assert len(self.records) == len(set([record.pos for record in self.records])), "There must 1 and only 1 record per position! "

        self.__find_calls()

        self.__get_variants()

        self.snp_distance = numpy.sum(self.is_snp)

    def __repr__(self):
        '''Overload the print function to write a summary of the VCF file

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
        '''
        Private method to find changes within the genome based on the variant file.

        Creates calls dict used elsewhere.
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

            # if the REF, ALT pair are the same length, check if we can decompose into SNPs
            if len(record.ref)==len(variant):

                for counter,(before,after) in enumerate(zip(record.ref,variant)):

                    # only make a change if the ALT  base is different to the REF base
                    if before!=after:

                        metadata={}
                        metadata['call']=after
                        metadata['ref']=before
                        metadata['pos']=counter
                        vcf_info={}
                        vcf_info=copy.deepcopy(record.values)
                        vcf_info['REF']=record.ref
                        vcf_info['ALTS']=record.alts
                        metadata['original_vcf_row']=vcf_info

                        self.calls[(index+counter,variant_type)]=metadata

            # otherwise the REF, ALT pair are different lengths
            else:

                # try and simplify into a single indel and multiple SNPs
                mutations = self._simplify_call(record.ref, variant)

                for (p, type_, bases) in mutations:

                    # p = max(p - 1, 0)
                    if type_ in ["ins", "del"]:

                        indel_length = len(bases)
                        # if type_ == "del":
                        #     indel_length *= -1
                        metadata = {}
                        metadata['call'] = (type_,bases)
                        metadata['ref']=record.ref[0]
                        metadata['pos']=p
                        vcf_info={}
                        vcf_info=copy.deepcopy(record.values)
                        vcf_info['REF']=record.ref
                        vcf_info['ALTS']=record.alts
                        metadata['original_vcf_row']=vcf_info

                        self.calls[(index+p,'indel')]=metadata

                    else:

                        metadata = {}
                        metadata['call'] = bases[1]
                        metadata['ref'] = bases[0]
                        metadata['pos'] = p
                        vcf_info = {}
                        vcf_info = copy.deepcopy(record.values)
                        vcf_info['REF'] = record.ref
                        vcf_info['ALTS'] = record.alts
                        metadata['original_vcf_row'] = vcf_info

                        self.calls[(index+p,variant_type)] = metadata

    def _simplify_call(self, ref, alt):
        '''Private method to simplify a complex call into one indel and multiple SNPs.

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
        if indel == "ins":
            #Ins after index, del at index so adjust ins
            start -= 1
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
        '''Convert the VCFFile to a pandas DataFrame.

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
        filter_ = []
        values = {}
        for record in self.records:
            chroms.append(record.chrom)
            pos.append(record.pos)
            refs.append(record.ref)
            alts.append(record.alts)
            qual.append(record.qual)
            infos.append(record.info)
            filter_.append(record.filter)
            for key in record.values:
                if values.get(key) is None:
                    values[key] = [record.values[key]]
                else:
                    values[key].append(record.values[key])
        df = pandas.DataFrame({
            "CHROM": chroms, "POS": pos,
            "REF": refs, "ALTS": alts,
            "QUAL": qual, "INFO": infos, "FILTER": filter_,
            **values
            })
        df.attrs = meta_data
        return df

    def __get_variants(self):
        '''
        Private method to pull the variants out of the VCFFile object.

        Builds arrays of the variant calls, and their respective genome indices, as well as
        masks to show whether there is a snp, het, null or indel call at the corresponding genome index:
        i.e is_snp[genome.nucleotide_number == indices[i]] gives a bool to determine if a genome has a SNP call at this position
        '''

        alts=[]
        variants = []
        indices = []
        refs=[]
        positions=[]
        is_snp = []
        is_het = []
        is_null = []
        is_indel = []
        indel_length = []
        metadata = defaultdict(list)

        for (index,type_) in sorted(list(self.calls.keys())):
            indices.append(index)
            call = self.calls[(index,type_)]['call']
            alt=call
            ref = self.calls[(index,type_)]['ref']
            if hasattr(self, 'genome'):
                assert self.genome.nucleotide_sequence[self.genome.nucleotide_index==index]==ref, 'reference nucleotide in VCF does not match the supplied genome at index position '+str(index)
            refs.append(ref)
            pos = self.calls[(index,type_)]['pos']
            positions.append(pos)
            #Update the masks with the appropriate types
            if type_ == 'indel':
                #Convert to ins_x or del_x rather than tuple
                variant = str(index)+"_"+call[0]+"_"+str(call[1])
                alt=call[1]
                is_indel.append(True)
                if call[0]=='ins':
                    indel_length.append(len(call[1]))
                else:
                    indel_length.append(-1*len(call[1]))
                is_snp.append(False)
                is_het.append(False)
                is_null.append(False)
            elif type_ == "snp":
                variant = str(index)+ref+'>'+call
                is_indel.append(False)
                indel_length.append(0)
                is_snp.append(True)
                is_het.append(False)
                is_null.append(False)
            elif type_ == 'het':
                variant = str(index)+ref+'>'+alt
                is_indel.append(False)
                indel_length.append(0)
                is_snp.append(False)
                is_het.append(True)
                is_null.append(False)
            elif type_ == 'null':
                variant = str(index)+ref+'>'+alt
                is_indel.append(False)
                indel_length.append(0)
                is_snp.append(False)
                is_het.append(False)
                is_null.append(True)
            alts.append(alt)
            variants.append(variant)
            for key in self.calls[(index,type_)]['original_vcf_row']:
                metadata[key].append(self.calls[(index,type_)]['original_vcf_row'][key])

        #Convert to numpy arrays for neat indexing
        self.alt_nucleotides=numpy.array(alts)
        self.variants = numpy.array(variants)
        self.nucleotide_index = numpy.array(indices)
        self.is_indel = numpy.array(is_indel)
        self.indel_length=numpy.array(indel_length)
        self.is_snp = numpy.array(is_snp)
        self.is_het = numpy.array(is_het)
        self.is_null = numpy.array(is_null)
        self.ref_nucleotides=numpy.array(refs)
        self.pos=numpy.array(positions)
        self.metadata = dict()
        for key in metadata:
            self.metadata[key] = numpy.array(metadata[key], dtype=object)
