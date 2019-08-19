
import numpy

# FIXME: problems with rrs, mfpB

class Gene(object):

    """Gene object that uses underlying numpy arrays"""

    def __init__(self,gene_name=None,sequence=None,index=None,position=None,codes_protein=False,first_amino_acid_position=1,first_nucleotide_index=None,promoter_sequence=None):

        assert gene_name is not None, "must provide a gene name!"
        self.gene_name=gene_name

        assert codes_protein in [True,False], gene_name+": coding_region must be True or False!"

        assert isinstance(sequence,numpy.ndarray), gene_name+": open reading frame sequence must be a Numpy array!"
        assert isinstance(index,numpy.ndarray), gene_name+": genome indices must be a Numpy array!"
        assert isinstance(position,numpy.ndarray), gene_name+": gene positions must be a Numpy array!"

        assert numpy.issubdtype(position.dtype.type,numpy.integer)
        assert numpy.issubdtype(index.dtype.type,numpy.integer)

        sequence=numpy.char.lower(sequence)
        assert numpy.count_nonzero(numpy.isin(sequence,['a','t','c','g','x','z']))==len(sequence)

        promoter_mask=position<0
        cds_mask=position>0

        if position[0]<position[-1]:

            self.reverse_strand=False
            self.gene_sequence=sequence[cds_mask]
            self.gene_position=position[cds_mask]
            self.gene_index=index[cds_mask]
            self.promoter_sequence=sequence[promoter_mask]
            self.promoter_position=position[promoter_mask]
            self.promoter_index=index[promoter_mask]
        else:
            self.reverse_strand=True
            self.gene_sequence=sequence[cds_mask][::-1]
            self.gene_position=position[cds_mask][::-1]
            self.gene_index=index[cds_mask][::-1]
            self.promoter_sequence=sequence[promoter_mask][::-1]
            self.promoter_position=position[promoter_mask][::-1]
            self.promoter_index=index[promoter_mask][::-1]

        self.gene_number_nucleotides=len(self.gene_sequence)
        self.promoter_number_nucleotides=len(self.promoter_sequence)
        self.total_number_nucleotides=len(sequence)

        # self.orf_first_nucleotide_index=self.orf_index[0]-1
        # self.orf_last_nucleotide_index=self.orf_index[-1]

        self.codes_protein=codes_protein

        self._setup_conversion_dicts()

        if self.codes_protein:
            # assert (self.gene_number_nucleotides%3)==0, gene_name+": must have an exact multiple of three bases in the open reading frame sequence array!"
            # remaining_bases=self.gene_number_nucleotides%3

            # if remaining_bases>0:
            #     self.gene_sequence=self.gene_sequence[:-remaining_bases]
            #     self.gene_position=self.gene_position[:-remaining_bases]
            #     self.gene_index=self.gene_index[:-remaining_bases]
            self._translate_sequence()
        else:
            self.amino_acids=None
            self.amino_acids_position=None
            self.triplets=None

    def _translate_sequence(self):

        self.amino_acids_position=numpy.unique(self.gene_position)

        tmp=[]
        for resid in self.amino_acids_position:
            triplet=''.join(i for i in self.gene_sequence[self.gene_position==resid])
            if len(triplet)==3:
                tmp.append(triplet)

        self.triplets=numpy.array(tmp)

        # now translate the triplets into amino acids using this new dictionary
        self.amino_acids=numpy.array([self.triplet_to_amino_acid[i] for i in self.triplets])

    def _setup_conversion_dicts(self):

        bases = ['t', 'c', 'a', 'g', 'x', 'z']
        aminoacids = 'FFLLXZSSSSXZYY!!XZCC!WXZXXXXXXZZZZXZLLLLXZPPPPXZHHQQXZRRRRXZXXXXXXZZZZXZIIIMXZTTTTXZNNKKXZSSRRXZXXXXXXZZZZXZVVVVXZAAAAXZDDEEXZGGGGXZXXXXXXZZZZXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZZZZZXZZZZZXZZZZZXZXXXXXXZZZZXZ'
        self.codons = numpy.array([a+b+c for a in bases for b in bases for c in bases])
        self.triplet_to_amino_acid = dict(zip(self.codons, aminoacids))
        self.amino_acids_of_codons=numpy.array([self.triplet_to_amino_acid[i] for i in self.codons])

    def __repr__(self):

        string_length=5

        output=self.gene_name+" gene\n"
        output+="%i nucleotides" % self.total_number_nucleotides
        if self.codes_protein:
            output+=", codes for protein\n"
        else:
            output+="\n"
        if self.promoter_sequence.size!=0:
            output+="".join(i for i in self.promoter_sequence[:string_length])
            output+="..."
            output+="".join(i for i in self.promoter_sequence[-string_length:])
            output+="\n"
            output+="".join(str(i)+" " for i in self.promoter_position[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in self.promoter_position[-string_length:])
            output+="\n"
        else:
            output+="promoter found in adjacent gene(s)\n"
        if self.codes_protein:
            output+="".join(i for i in self.amino_acids[:string_length])
            output+="..."
            output+="".join(i for i in self.amino_acids[-string_length:])
            output+="\n"
            output+="".join(str(i)+" " for i in self.amino_acids_position[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in self.amino_acids_position[-string_length:])
            output+="\n"
        else:
            output+="".join(i for i in self.gene_sequence[:string_length])
            output+="..."
            output+="".join(i for i in self.gene_sequence[-string_length:])
            output+="\n"
            output+="".join(str(i)+" " for i in self.gene_position[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in self.gene_position[-string_length:])
            output+="\n"

        if output.strip()=="":
            output=None

        return(output)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genomes
        """

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.gene_name==other.gene_name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        promoter_mask=self.promoter_sequence!=other.promoter_sequence
        if self.codes_protein:
            gene_mask=self.amino_acids!=other.amino_acids
            tmp=self.amino_acids_position[gene_mask]
        else:
            gene_mask=self.gene_sequence!=other.gene_sequence
            tmp=self.gene_position[gene_mask]

        return(self.promoter_index[promoter_mask],tmp)

    def valid_element(self, element=None):

        return(True)
