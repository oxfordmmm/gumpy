import pysam


class VCFRecord(object):
    '''
    Class for VCF records
    '''
    def __init__(self, record, sample):
        '''Constructor, pulls the data out of the pysam object
        into a more useable format

        Args:
            record (pysam.libcbcf.VariantRecord): The record object
            sample_num (str) : Name of the sample to consider. Used for possible cases
                                where there is more than 1 sample per record
        '''
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
    '''
    def __init__(self, filename):
        '''
        Constructor. Reads the VCF file and sets up values
        Args:
            filename (str) : The name of the VCF file
        '''
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

        self.find_changes()
    
    def find_changes(self):
        '''Function to find changes within the genome based on the variant file
        '''
        #Use a dict to store the changes with keys as indicies
        #{index: (base, {calls: base}) ...} where a base of * represents wildcard
        self.changes = {}
        for record in self.records:
            #VCF files are 1 indexed so alter index to be 0 indexed
            key = record.pos - 1
            ref = record.ref
            alts = record.alts
            if len(alts) == 1:
                if len(alts[0]) > 1:
                    #Indel so use x
                    call = 'x'
                else:
                    #Single call so set it as such
                    call = alts[0]
            else:
                #Het call so use z
                call = 'z'

            #Get the number of reads for each call (and wildcard)
            n_reads = record.values.get("COV")
            #Ensure that there is a number of reads set
            assert n_reads is not None, "There is no number of reads on at least one line of the file (COV)"

            #Get the most frequent call
            if max(n_reads) == 0:
                #Wildcard is the most common (TODO)
                # assert False, "Wildcard is the most common..."
                call = 'x'
            alts = ('*', ) + alts
            self.changes[key] = (call, {n: base for (n, base) in zip(n_reads, alts)})
        # print(self.changes)










'''
class VariantFile(object):

    def __init__(self, filename=None):
        a=2


class VCFFile(VariantFile):

    def __init__(self, filepath):

        self.filepath=pathlib.Path(filepath)
        # read and store the header information from the VCF file
        self._read_header()


    def _read_header(self):

        vcf_reader = pysam.VariantFile(self.filepath)

        # store the VCF version
        self.version=vcf_reader.header.version

        # store the reference contigs, allowing for multiple entries
        self.contigs={}

        for i in list(vcf_reader.header.contigs):
            self.contigs[i]={ 'length': int(vcf_reader.header.contigs[i].length)}

        # store the filters applied, allowing for multiple filters
        self.filters={}
        for i in list(vcf_reader.header.filters):
            self.filters[i]={ 'description': vcf_reader.header.filters[i].description}

        self.info_fields={}

        for i in list(vcf_reader.header.formats):

            self.info_fields[i] = {  'type':vcf_reader.header.formats[i].type.lower(),\
                                    'description':vcf_reader.header.formats[i].description,\
                                    'number':vcf_reader.header.formats[i].number }

    def __repr__(self):
        output=str(self.filepath)+'\n'
        for i in self.contigs:
            output+="%s, length %i" % (i,self.contigs[i]['length'])            
        return(output)

    def apply_to_genome(self, genome):

        assert type(genome)==gumpy.Genome, 'the passed object must be an instance of a gumpy.Genome object'

        assert len(self.contigs)==1, 'only VCF files that apply to a single contig can be parsed'

        for ref_name in self.contigs:
            assert genome.genome_length==self.contigs[ref_name]['length'], 'recorded length of reference does not match passed Genome!'

        vcf_reader = pysam.VariantFile(self.filepath)

        counter=0

        for record in vcf_reader:

            for sample_idx, (sample_name, sample_info) in enumerate(
                record.samples.items()
            ):

                print(sample_name, sample_info.keys())

                counter+=1

            if counter>10:
                break
'''