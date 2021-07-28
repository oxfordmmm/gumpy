import pysam, numpy
import pandas as pd
from gumpy import VCFDifference


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
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
            return
        else:
            record = args[0]
            sample = args[1]
        #Save some of the easier to get to attributes
        self.chrom = record.chrom
        self.pos = record.pos
        self.ref = record.ref
        self.alts = record.alts
        self.qual = record.qual

        #Get the filter attribute value
        assert len(record.filter.items()) >= 0, "A record has more than 1 filter set!"
        if len(record.filter.items()) == 0:
            self.filter = "."
        else:
            self.filter = record.filter.items()[0][0]
        
        #Get the info field
        self.info = {}
        for (key, value) in record.info.items():
            self.info[key] = value

        #Get the values
        self.values = {}
        for (key, item) in record.samples[sample].items():
            if key == "GT_CONF":
                #Due to how pysam reads floats, there are some erroneously long dps, so round
                item = round(item, 2)
            self.values[key] = item
    
    def __repr__(self):
        '''Pretty print the record
        '''
        s = self.chrom + "\t"
        s += str(self.pos) + "\t"
        s += self.ref + "\t"
        s += str(self.alts) + "\t"
        s += str(self.qual) + "\t"
        s += self.filter + "\t"
        for val in self.values.keys():
            s += str(val)+":"
        s = s[:-1] + "\t"
        for val in self.values.values():
            s += str(val)+":"
        s = s[:-1]
        s += "\n"
        return s

class VariantFile(object):
    '''
    Class to instanciate a variant file (VCF)
    Used to apply a VCF file to a genome
    Instance variables:
        `VCF_VERSION` (tuple(int)): Tuple of ints to show the VCF version of the file. e.g 3.2 would be (3, 2).
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
            allowed_kwargs = ['VCF_VERSION', 'contig_lengths', 'formats', 'records', 'changes']
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
            return
        else:
            filename = args[0]
        #Use pysam to parse the VCF
        vcf = pysam.VariantFile(filename)

        #Get some basic metadata
        self.VCF_VERSION = vcf.version
        
        #Get the contig lengths from the header
        self.contig_lengths = {}
        for name in list(vcf.header.contigs):
            self.contig_lengths[name] = vcf.header.contigs[name].length
        # print(self.contig_lengths)

        #Get the formats
        self.formats = {}
        for format in vcf.header.formats.keys():
            description = vcf.header.formats[format].description
            id = vcf.header.formats[format].id
            f_type = vcf.header.formats[format].type
            self.formats[format] = {
                "description" : description,
                "id": id,
                "type": f_type
            }
        # print(self.formats)

        #Get the records
        self.records = []
        for record in list(vcf):
            for sample in record.samples.keys():
                self.records.append(VCFRecord(record, sample))
        
        # print(self.records)

        self.__find_changes()
    
    def __find_changes(self):
        '''Function to find changes within the genome based on the variant file
        '''
        #Use a dict to store the changes with keys as indicies
        #{index: (base, {calls: base}) ...} where a base of * represents wildtype
        #Where index is 0 indexed for array access
        self.changes = {}
        # i = 0
        for record in self.records:
            #VCF files are 1 indexed so alter index to be 0 indexed
            key = record.pos - 1
            ref = record.ref
            alts = record.alts
            #Check for singular null call
            if alts is None:
                alts = ('x', )
            alts = [call if call is not None else 'x' for call in alts]
            alts = tuple([call.lower() if len(call) == 1 else numpy.array(list(call.lower())) for call in alts])
            if len(alts) == 1:
                if len(alts[0]) > 1:
                    #Indel so use x
                    call = alts
                else:
                    #Single call so set it as such
                    call = alts[0]
            else:
                '''
                4 Different types of het calls found:
                    1. SNP het call e.g A->(G, T)
                    2. indel het call e.g A->(ACT, AGT)
                    3. indel-indel het call e.g ACG->(CGT, ATG)
                    4. Mixed SNP and indel het call e.g A->(A, T, ACT, GTC)
                How to determine insertion/deletion from an indel??
                '''
                if numpy.any([len(call) > 1 for call in alts]):
                    #There is at least one indel call so record that
                    # call = ("@@", ) +alts
                    call = 'z'
                elif len(set(alts)) == len(alts) or len(set(alts)) > 1:
                    #Check to make sure the multi-calling has different calls
                    #Het call so use z
                    call = 'z'
                else:
                    assert False, "Weird call: "+str(alts)
                

            #Get the number of reads for each call (and wildcard)
            n_reads = record.values.get("COV")
            #Ensure that there is a number of reads set
            assert n_reads is not None, "There is no number of reads on at least one line of the file (COV)"

            #Get the most frequent call
            # if max(n_reads) == 0:
            #     #Wildtype is the most common (TODO)
            #     call = 'x'
            alts = ('*', ) + alts
            self.changes[key] = (call, list(zip(n_reads, alts)))

    def to_df(self):
        '''Convert the VCFRecord to a pandas DataFrame. 
        Metadata is stored in the `attrs` attribute of the DataFrame which may break with some operations 
        (but pandas does not currently have a robust method for metadata storage...)

        Returns:
            pandas.DataFrame: DataFrame containing all of the information from the VCF file
        '''        
        meta_data = {
            "VCF_VERSION": self.VCF_VERSION,
            "contig_lengths": self.contig_lengths,
            "formats": self.formats
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
    
    def difference(self, genome):
        '''Takes in a Genome object, returning an object detailing the full differences caused by this VCF including amino acid differences

        Args:
            genome (gumpy.Genome): A reference Genome object
        '''        
        return VCFDifference(self, genome)