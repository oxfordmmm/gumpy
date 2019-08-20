
import numpy

# FIXME: problems with rrs, mfpB

class Gene(object):

    """Gene object that uses underlying numpy arrays"""

    def __init__(self,gene_name=None,sequence=None,index=None,numbering=None,codes_protein=True,feature_type=None):

        assert gene_name is not None, "must provide a gene name!"
        self.gene_name=gene_name

        assert codes_protein in [True,False], gene_name+": codes_protein must be True or False!"
        self.codes_protein=codes_protein

        assert isinstance(sequence,numpy.ndarray), gene_name+": sequence of bases must be a Numpy array!"

        assert isinstance(index,numpy.ndarray), gene_name+": genome indices must be a Numpy array of integers!"
        assert numpy.issubdtype(index.dtype.type,numpy.integer), gene_name+": genome indices must be a Numpy array of integers!"

        assert isinstance(numbering,numpy.ndarray), gene_name+": gene numbering must be a Numpy array of integers!"
        assert numpy.issubdtype(numbering.dtype.type,numpy.integer), gene_name+": gene numbering must be a Numpy array of integers!"

        sequence=numpy.char.lower(sequence)
        assert numpy.count_nonzero(numpy.isin(sequence,['a','t','c','g','x','z']))==len(sequence), gene_name+": sequence can only contain a,t,c,g,z,x"

        promoter_mask=numbering<0
        cds_mask=numbering>0

        if numbering[0]<numbering[-1]:
            self.on_noncoding_strand=False
            self.sequence=sequence
            self.numbering=numbering
            self.index=index
            self.is_cds=cds_mask
            self.is_promoter=promoter_mask
        else:
            self.on_noncoding_strand=True
            self.sequence=sequence[::-1]
            self.numbering=numbering[::-1]
            self.index=index[::-1]
            self.is_cds=cds_mask[::-1]
            self.is_promoter=promoter_mask[::-1]

        self.cds_index_start=numpy.min(self.index[self.is_cds])
        self.cds_index_end=numpy.max(self.index[self.is_cds])
        self.cds_number_nucleotides=len(self.sequence[self.is_cds])

        if numpy.count_nonzero(self.is_promoter)>0:
            self.promoter_index_start=numpy.min(self.index[self.is_promoter])
            self.promoter_index_end=numpy.max(self.index[self.is_promoter])
            self.promoter_number_nucleotides=len(self.sequence[self.is_promoter])
        else:
            self.promoter_index_start=None
            self.promoter_index_end=None
            self.promoter_number_nucleotides=0

        self.total_number_nucleotides=len(sequence)


        if self.codes_protein:

            self._setup_conversion_dicts()
            self._translate_sequence()

        else:

            self.amino_acid_sequence=None
            self.amino_acid_numbering=None
            self.codons=None

    def _translate_sequence(self):

        # this will ensure that only amino acids with all three bases present
        unique,counts=numpy.unique(self.numbering[self.is_cds],return_counts=True)
        self.amino_acid_numbering=unique[counts==3]

        # try to optimize!
        shorter_numbering=self.numbering[self.is_cds]
        shorter_sequence=self.sequence[self.is_cds]
        trip=[]
        for resid in self.amino_acid_numbering:
            triplet=''.join(i for i in shorter_sequence[shorter_numbering==resid])
            trip.append(triplet)

        self.codons=numpy.array(trip)

        # now translate the triplets into amino acids using this new dictionary
        self.amino_acid_sequence=numpy.array([self.codon_to_amino_acid[i] for i in self.codons])

    def _setup_conversion_dicts(self):

        bases = ['t', 'c', 'a', 'g', 'x', 'z']
        aminoacids = 'FFLLXZSSSSXZYY!!XZCC!WXZXXXXXXZZZZXZLLLLXZPPPPXZHHQQXZRRRRXZXXXXXXZZZZXZIIIMXZTTTTXZNNKKXZSSRRXZXXXXXXZZZZXZVVVVXZAAAAXZDDEEXZGGGGXZXXXXXXZZZZXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZZZZZXZZZZZXZZZZZXZXXXXXXZZZZXZ'
        all_codons = numpy.array([a+b+c for a in bases for b in bases for c in bases])
        self.codon_to_amino_acid = dict(zip(all_codons, aminoacids))
        # self.amino_acids_of_codons=numpy.array([self.codon_to_amino_acid[i] for i in all_codons])

    def __repr__(self):

        string_length=5

        output=self.gene_name+" gene\n"
        output+="%i nucleotides" % self.total_number_nucleotides
        if self.codes_protein:
            output+=", codes for protein\n"
        else:
            output+="\n"
        promoter_sequence=self.sequence[self.is_promoter]
        if promoter_sequence.size!=0:
            output+="".join(i for i in promoter_sequence[:string_length])
            output+="..."
            output+="".join(i for i in promoter_sequence[-string_length:])
            output+="\n"
            promoter_numbering=self.numbering[self.is_promoter]
            output+="".join(str(i)+" " for i in promoter_numbering[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in promoter_numbering[-string_length:])
            output+="\n"
        else:
            output+="promoter likely in adjacent gene(s)\n"
        if self.codes_protein:
            output+="".join(i for i in self.amino_acid_sequence[:string_length])
            output+="..."
            output+="".join(i for i in self.amino_acid_sequence[-string_length:])
            output+="\n"
            output+="".join(str(i)+" " for i in self.amino_acid_numbering[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in self.amino_acid_numbering[-string_length:])
            output+="\n"
        else:
            gene_sequence=self.sequence[self.is_cds]
            output+="".join(i for i in gene_sequence[:string_length])
            output+="..."
            output+="".join(i for i in gene_sequence[-string_length:])
            output+="\n"
            gene_numbering=self.numbering[self.is_cds]
            output+="".join(str(i)+" " for i in gene_numbering[:string_length])
            output+="..."
            output+="".join(str(i)+" " for i in gene_numbering[-string_length:])
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

        mutations=[]
        if self.codes_protein:
            mask=self.amino_acid_sequence!=other.amino_acid_sequence
            pos=self.amino_acid_numbering[mask]
            ref=self.amino_acid_sequence[mask]
            alt=other.amino_acid_sequence[mask]

            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(p)+a)

            mask=(self.sequence!=other.sequence) & self.is_promoter
            pos=self.numbering[mask]
            ref=self.sequence[mask]
            alt=other.sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(p)+a)
        else:
            mask=self.sequence!=other.sequence
            pos=self.numbering[mask]
            ref=self.sequence[mask]
            alt=other.sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(p)+a)

        if not mutations:
            mutations=None
        return(mutations)

    def valid_element(self, element=None):

        return(True)
