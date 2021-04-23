

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

            if record.type in ['CDS','rRNA']:

                gene_name=None
                type=None
                codes_protein=True

                if 'gene' in record.qualifiers.keys():

                    gene_name=record.qualifiers['gene'][0]

                    if record.type=='rRNA':
                        type="RNA"
                    else:
                        type="GENE"
                        codes_protein=True

                elif 'locus_tag' in record.qualifiers.keys():

                    gene_name=record.qualifiers['locus_tag'][0]

                    if record.type=='rRNA':
                        type="RNA"
                    else:
                        type="LOCUS"
                        codes_protein=True

                else:
                    continue

                if gene_name is not None and (gene_subset is None or gene_name in gene_subset):

                    assert len(gene_name)<=20, "Gene "+gene_name+" is too long at "+str(len(gene_name))+" chars; need to change numpy.zeros definiton U20 to match"
