

import numpy

from Bio import SeqIO

class NewGenome(object):

    def __init__(self,\
                genbank_file,\
                show_progress_bar=False,\
                gene_subset=None,\
                default_promoter_length=100):

        assert default_promoter_length>=0, "the promoter length must be a positive integer!"
        assert isinstance(default_promoter_length,int), "the promoter length must be a positive integer!"

        reference_genome=SeqIO.read(genbank_file,'genbank')

        # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
        self.genome_sequence=numpy.array([i.lower() for i in str(reference_genome.seq)])

        self.name=reference_genome.name
        self.id=reference_genome.id
        self.description=reference_genome.description

        # store the length of the genome
        self.genome_length=len(self.genome_sequence)

        # create an array of the genome indices
        self.genome_index=numpy.arange(1,self.genome_length+1,dtype="int")

        assert len(reference_genome.annotations['accessions'])==1, 'only GenBank files with a single accessions currently allowed'

        self.annotations={}
        for i in reference_genome.annotations.keys():
            self.annotations[i]=reference_genome.annotations[i]

        for record in reference_genome.features:

            if record.type not in ['CDS','rRNA']:
                break

            gene_name=None
            type=None
            codes_protein=True

            if 'gene' in record.qualifiers.keys():
                gene_name=record.qualifiers['gene'][0]
                type='GENE'

            elif 'locus_tag' in record.qualifiers.keys():
                gene_name=record.qualifiers['locus_tag'][0]
                type="LOCUS"

            if record.type=='rRNA':
                type="RNA"
                codes_protein=False

            if gene_name is None or (gene_subset is not None or gene_name not in gene_subset):
                break

            assert len(gene_name)<=20, "Gene "+gene_name+" is too long at "+str(len(gene_name))+" chars; need to change numpy.zeros definiton U20 to match"

            gene_start=int(record.location.start)
            gene_end=int(record.location.end)

            gene_name+="_2" if gene_name in genes_found_so_far else ''
            genes_found_so_far.append(gene_name)

            if record.strand==1:

                if previous_gene_reversed:
                    gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end) & (~self.genome_is_cds)
                else:
                    gene_mask=(self.genome_index>gene_start) & (self.genome_index<=gene_end)

                if codes_protein:
                    gene_coding_numbering=numpy.floor_divide(self.genome_index[gene_mask]-gene_start+2,3)
                else:
                    gene_coding_numbering=self.genome_index[gene_mask]-gene_start

                gene_coding_position=self.genome_index[gene_mask]-gene_start

                previous_gene_reversed=False


            self.genome_feature_type[mask]=type
            self.genome_feature_name[mask]=gene_name

            self.genome_position[gene_mask]=gene_coding_numbering

            if codes_protein:
                self.genome_amino_acid_number[gene_mask]=gene_coding_numbering

            self.genome_nucleotide_number[gene_mask]=gene_coding_position

            self.genome_is_cds[gene_mask]=True

            self._gene_type[gene_name]=type
            self._gene_codes_protein[gene_name]=codes_protein
