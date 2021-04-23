import pathlib, pysam, gumpy

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
        
