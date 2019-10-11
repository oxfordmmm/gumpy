
import numpy, pandas

# FIXME: problems with rrs, mfpB

class Gene(object):

    """Gene object that uses underlying numpy arrays"""

    def __init__(self,gene_name=None,sequence=None,index=None,numbering=None,positions=None,is_indel=None,indel_length=None,codes_protein=True,feature_type=None):

        assert gene_name is not None, "must provide a gene name!"
        self.gene_name=gene_name

        self.gene_type=feature_type
        assert codes_protein in [True,False], gene_name+": codes_protein must be True or False!"
        self.codes_protein=codes_protein

        assert isinstance(sequence,numpy.ndarray), gene_name+": sequence of bases must be a Numpy array!"

        assert isinstance(index,numpy.ndarray), gene_name+": genome indices must be a Numpy array of integers!"
        assert numpy.issubdtype(index.dtype.type,numpy.integer), gene_name+": genome indices must be a Numpy array of integers!"

        assert isinstance(numbering,numpy.ndarray), gene_name+": gene numbering must be a Numpy array of integers!"
        # assert numpy.issubdtype(numbering.dtype.type,numpy.float64), gene_name+": gene numbering must be a Numpy array of integers!"

        sequence=numpy.char.lower(sequence)
        assert numpy.count_nonzero(numpy.isin(sequence,['a','t','c','g','x','z']))==len(sequence), gene_name+": sequence can only contain a,t,c,g,z,x"

        promoter_mask=numbering<0
        cds_mask=numbering>0

        if numbering[0]<numbering[-1]:
            self.on_noncoding_strand=False
            self.sequence=sequence
            self.numbering=numbering
            self.positions=positions
            self.is_indel=is_indel
            self.indel_length=indel_length
            self.index=index
            self.is_cds=cds_mask
            self.is_promoter=promoter_mask
        else:
            self.on_noncoding_strand=True
            self.sequence=sequence[::-1]
            self.numbering=numbering[::-1]
            self.positions=positions[::-1]
            self.is_indel=is_indel[::-1]
            self.indel_length=indel_length[::-1]
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
            self.amino_acid_numbering=0
            self.codons=None

    def _translate_sequence(self):

        # this will ensure that only amino acids with all three bases present
        unique,counts=numpy.unique(self.numbering[self.is_cds],return_counts=True)
        self.amino_acid_numbering=unique[counts==3]
        self.amino_acid_numbering=self.amino_acid_numbering.astype(int)

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

    def list_mutations_wrt(self,other):

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.gene_name==other.gene_name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        mutations=[]
        if self.codes_protein:

            # deal first with the coding sequence, which we will view at the amino acid level
            # using codons rather than amino acids will detect synonymous mutations as well
            mask=(self.codons!=other.codons)
            pos=self.amino_acid_numbering[mask]
            ref=other.amino_acid_sequence[mask]
            alt=self.amino_acid_sequence[mask]

            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)

            mask=(self.sequence!=other.sequence) & self.is_promoter
            pos=self.numbering[mask]
            ref=other.sequence[mask]
            alt=self.sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)
        else:
            mask=self.sequence!=other.sequence
            pos=self.positions[mask]
            ref=other.sequence[mask]
            alt=self.sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)

        mask=self.is_indel
        pos=list(self.positions[mask])
        length=list(self.indel_length[mask])
        for (p,l) in zip(pos,length):
            mutations.append(str(int(p))+"_indel")

        if not mutations:
            mutations=None

        return(mutations)


    def table_mutations_wrt(self,other):

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.gene_name==other.gene_name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        MUTATIONS_dict={}
        MUTATIONS_columns=['GENE','MUTATION','REF','ALT','POSITION','AMINO_ACID_NUMBER','GENOME_INDEX','NUCLEOTIDE_NUMBER','IS_SNP','IS_INDEL','IN_CDS','IN_PROMOTER','ELEMENT_TYPE','MUTATION_TYPE','INDEL_LENGTH','INDEL_1','INDEL_2']
        for cols in MUTATIONS_columns:
            MUTATIONS_dict[cols]=[]

        if self.codes_protein:

            # deal first with the coding sequence, which we will view at the amino acid level
            # using codons rather than amino acids will detect synonymous mutations as well
            mask=(self.codons!=other.codons)
            pos=self.amino_acid_numbering[mask]
            ref=other.amino_acid_sequence[mask]
            alt=self.amino_acid_sequence[mask]
            ref_codon=other.codons[mask]
            alt_codon=self.codons[mask]

            for (r,p,a,rc,ac) in zip(ref,pos,alt,ref_codon,alt_codon):
                mut=r+str(int(p))+a
                MUTATIONS_dict['GENE'].append(self.gene_name)
                MUTATIONS_dict['MUTATION'].append(mut)
                MUTATIONS_dict['REF'].append(rc)
                MUTATIONS_dict['ALT'].append(ac)
                MUTATIONS_dict['POSITION'].append(p)
                MUTATIONS_dict['AMINO_ACID_NUMBER'].append(p)
                MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(0)
                MUTATIONS_dict['GENOME_INDEX'].append(0)
                MUTATIONS_dict['IS_SNP'].append(True)
                MUTATIONS_dict['IS_INDEL'].append(False)
                MUTATIONS_dict['IN_CDS'].append(True)
                MUTATIONS_dict['IN_PROMOTER'].append(False)
                MUTATIONS_dict['INDEL_LENGTH'].append(0)
                MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
                MUTATIONS_dict['MUTATION_TYPE'].append('AAM')
                MUTATIONS_dict['INDEL_1'].append(None)
                MUTATIONS_dict['INDEL_2'].append(None)

            mask=(self.sequence!=other.sequence) & self.is_promoter
            pos=self.numbering[mask]
            ref=other.sequence[mask]
            alt=self.sequence[mask]
            idx=self.index[mask]
            for (r,p,a,i) in zip(ref,pos,alt,idx):
                mut=r+str(int(p))+a
                MUTATIONS_dict['GENE'].append(self.gene_name)
                MUTATIONS_dict['MUTATION'].append(mut)
                MUTATIONS_dict['REF'].append(r)
                MUTATIONS_dict['ALT'].append(a)
                MUTATIONS_dict['POSITION'].append(p)
                MUTATIONS_dict['AMINO_ACID_NUMBER'].append(0)
                MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
                MUTATIONS_dict['GENOME_INDEX'].append(i)
                MUTATIONS_dict['IS_SNP'].append(True)
                MUTATIONS_dict['IS_INDEL'].append(False)
                MUTATIONS_dict['IN_CDS'].append(False)
                MUTATIONS_dict['IN_PROMOTER'].append(True)
                MUTATIONS_dict['INDEL_LENGTH'].append(0)
                MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
                MUTATIONS_dict['MUTATION_TYPE'].append('SNP')
                MUTATIONS_dict['INDEL_1'].append(None)
                MUTATIONS_dict['INDEL_2'].append(None)

        else:
            mask=self.sequence!=other.sequence
            pos=self.positions[mask]
            ref=other.sequence[mask]
            alt=self.sequence[mask]
            idx=self.index[mask]
            for (r,p,a,i) in zip(ref,pos,alt,idx):
                if p<0:
                    is_promoter=True
                    is_cds=False
                else:
                    is_promoter=False
                    is_cds=True
                mut=r+str(int(p))+a
                MUTATIONS_dict['GENE'].append(self.gene_name)
                MUTATIONS_dict['MUTATION'].append(mut)
                MUTATIONS_dict['REF'].append(r)
                MUTATIONS_dict['ALT'].append(a)
                MUTATIONS_dict['POSITION'].append(p)
                MUTATIONS_dict['AMINO_ACID_NUMBER'].append(0)
                MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
                MUTATIONS_dict['GENOME_INDEX'].append(i)
                MUTATIONS_dict['IS_SNP'].append(True)
                MUTATIONS_dict['IS_INDEL'].append(False)
                MUTATIONS_dict['IN_CDS'].append(is_cds)
                MUTATIONS_dict['IN_PROMOTER'].append(is_promoter)
                MUTATIONS_dict['INDEL_LENGTH'].append(0)
                MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
                MUTATIONS_dict['MUTATION_TYPE'].append('SNP')
                MUTATIONS_dict['INDEL_1'].append(None)
                MUTATIONS_dict['INDEL_2'].append(None)


        mask=self.is_indel
        pos=list(self.positions[mask])
        num=list(self.numbering[mask])
        length=list(self.indel_length[mask])
        idx=self.index[mask]
        for (p,l,n,i) in zip(pos,length,num,idx):
            if p<0:
                is_promoter=True
                is_cds=False
                n=None
            else:
                is_promoter=False
                is_cds=True
            mut0=str(int(p))+"_indel"
            if l>0:
                mut1=str(int(p))+"_ins"
                mut2=mut1+"_"+str(l)
            else:
                mut1=str(int(p))+"_del"
                mut2=mut1+"_"+str(-1*l)
            MUTATIONS_dict['GENE'].append(self.gene_name)
            MUTATIONS_dict['MUTATION'].append(mut0)
            MUTATIONS_dict['REF'].append(None)
            MUTATIONS_dict['ALT'].append(None)
            MUTATIONS_dict['POSITION'].append(p)
            MUTATIONS_dict['AMINO_ACID_NUMBER'].append(n)
            MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
            MUTATIONS_dict['GENOME_INDEX'].append(i)
            MUTATIONS_dict['IS_SNP'].append(False)
            MUTATIONS_dict['IS_INDEL'].append(True)
            MUTATIONS_dict['IN_CDS'].append(is_cds)
            MUTATIONS_dict['IN_PROMOTER'].append(is_promoter)
            MUTATIONS_dict['INDEL_LENGTH'].append(l)
            MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
            MUTATIONS_dict['MUTATION_TYPE'].append('INDEL')
            MUTATIONS_dict['INDEL_1'].append(mut1)
            MUTATIONS_dict['INDEL_2'].append(mut2)

        if len(MUTATIONS_dict['POSITION'])>0:

            MUTATIONS_table=pandas.DataFrame(data=MUTATIONS_dict)

            MUTATIONS_table=MUTATIONS_table[MUTATIONS_columns]

            MUTATIONS_table=MUTATIONS_table.astype({'GENE':'category',\
                                                    'MUTATION':'str',\
                                                    'REF':'str',\
                                                    'ALT':'str',\
                                                    'POSITION':'float',\
                                                    'AMINO_ACID_NUMBER':'float',\
                                                    'NUCLEOTIDE_NUMBER':'float',\
                                                    'GENOME_INDEX':'float',\
                                                    'IS_SNP':'bool',\
                                                    'IS_INDEL':'bool',\
                                                    'IN_CDS':'bool',\
                                                    'IN_PROMOTER':'bool',\
                                                    'INDEL_LENGTH':'float',\
                                                    'ELEMENT_TYPE':'category',\
                                                    'MUTATION_TYPE':'category',\
                                                    'INDEL_1':'str',\
                                                    'INDEL_2':'str'})

            # MUTATIONS_table=MUTATIONS_table.replace({   'POSITION':0,\
            #                                             'NUCLEOTIDE_NUMBER':0,\
            #                                             'AMINO_ACID_NUMBER':0,\
            #                                             'GENOME_INDEX':0,\
            #                                             'INDEL_LENGTH':0  }, numpy.nan)

            return(MUTATIONS_table)

        else:

            return(None)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genes
        """

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.gene_name==other.gene_name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        positions=[]
        # if self.codes_protein:
        #
        #     # using codons rather than amino acids will detect synonymous mutations as well
        #     mask=self.codons!=other.codons
        #     pos=list(self.positions[mask])
        #     positions=pos
        #
        #     mask=(self.sequence!=other.sequence) & self.is_promoter
        #     pos=list(self.positions[mask])
        #     if positions is not None:
        #         positions=positions+pos
        #     else:
        #         positions=pos
        # else:
        mask=self.sequence!=other.sequence
        positions=list(self.positions[mask])

        mask=self.is_indel
        pos=list(self.positions[mask])
        positions=positions+pos

        if not positions:
            return None

        return(numpy.array(positions))


    def valid_variant(self, variant):

        assert variant is not None, "variant must be specified! e.g. FIXME"

        # since the smallest variant is 'a1c'
        assert len(variant)>=3, "a variant must have at least 3 characters e.g. a3c"

        before=variant[0]
        after=variant[-1]

        assert before!=after, "before and after are identical hence this is not a variant!"

        try:
            position=int(variant[1:-1])
        except:
            raise TypeError("position "+variant[1:-1]+" is not an integer!")

        # check that the specified base is actually a base!
        assert before in ['c','t','g','a'], before+" is not a valid nucleotide!"

        # having checked that position is an integer, make a mask
        mask=self.positions==position

        # if the mask contains anything other than a single True, the position is not in the genome
        assert numpy.count_nonzero(mask)==1, "position specified not in the genome!"

        # confirm that the given base matches what is in the genome
        assert before==self.sequence[mask][0], "base in genome is "+self.sequence[mask][0]+" but specified base is "+before

        # check that the base to be mutated to is valid (z=het, ?=any other base according to the grammar)
        assert after in ['c','t','g','a','?','z'], after+" is not a valid nucleotide!"

        return True
