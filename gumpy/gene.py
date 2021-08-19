import numpy

# FIXME: problems with rrs, mfpB
class Gene(object):
    """Gene object that uses underlying numpy arrays"""
    def __init__(self, *args, **kwargs):
        '''Constructor for the Gene object.
        Args:
            name (str, optional): Name of the gene. Defaults to None
            nucleotide_sequence (numpy.array, optional): Numpy array of the nucleotide sequence. Defaults to None
            index (numpy.array, optional): Numpy array of the gene indices. Defaults to None
            nucleotide_number (numpy.array, optional): Numpy array of the gene numbering. Defaults to None
            is_cds (numpy.array, optional): Numpy array to act as a mask for whether given elements are codons. Defaults to None
            is_promoter (numpy.array, optional): Numpy array to act as a mask for whether given elements are promoters. Defaults to None
            is_indel (numpy.array, optional): Numpy array to act as a mask for whether given elements are indels. Defaults to None
            indel_length (numpy.array, optional): Numpy array denoting the lengths of the indels whenever is_indel==True. Defaults to None
            reverse_complement (boolean, optional): Boolean showing whether this gene is reverse complement. Defaults to False
            codes_protein (boolean, optional): Boolean showing whether this gene codes a protein. Defaults to True
            feature_type (str, optional): The name of the type of feature that this gene represents. Defaults to None
            ribosomal_shifts (list(int), optional): Indices of repeated bases due to ribosomal frame shifting. Defaults to []
        '''
        #Set the kwargs
        #UPDATE THIS AS REQUIRED
        allowed_kwargs = ['name', 'nucleotide_sequence', 'index', 'nucleotide_number', 'is_cds', 'is_promoter', 
                        'is_indel', 'indel_length', 'codes_protein', 'reverse_complement', 'feature_type', 
                        'triplet_number', 'total_number_nucleotides', 'codon_to_amino_acid', 'amino_acid_number', 
                        'codons', 'amino_acid_sequence', 'ribosomal_shifts']
        #If reloading a Gene, just set the attributes and return
        if "reloading" in kwargs.keys():
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
            return
        
        #Set default values based on kwargs
        name = kwargs.get("name")
        nucleotide_sequence = kwargs.get("nucleotide_sequence")
        index = kwargs.get("index")
        nucleotide_number = kwargs.get("nucleotide_number")
        is_cds = kwargs.get("is_cds")
        is_promoter = kwargs.get("is_promoter")
        is_indel = kwargs.get("is_indel")
        indel_length = kwargs.get("indel_length")
        reverse_complement = kwargs.get("reverse_complement", False)
        codes_protein = kwargs.get("codes_protein", True)
        feature_type = kwargs.get("feature_type")
        ribosomal_shifts = kwargs.get("ribosomal_shifts", [])

        assert name is not None, "must provide a gene name!"
        self.name=name

        self.feature_type=feature_type
        assert codes_protein in [True,False], name+": codes_protein must be True or False!"
        self.codes_protein=codes_protein

        assert reverse_complement in [True,False], name+": reverse_complement must be True or False!"
        self.reverse_complement=reverse_complement

        assert isinstance(nucleotide_sequence,numpy.ndarray), name+": sequence of bases must be a Numpy array!"

        assert isinstance(index,numpy.ndarray) and numpy.issubdtype(index.dtype.type,numpy.integer), name+": genome indices must be a Numpy array of integers!"

        assert isinstance(nucleotide_number,numpy.ndarray), name+": gene numbering must be a Numpy array of integers!"

        nucleotide_sequence=numpy.char.lower(nucleotide_sequence)
        assert numpy.count_nonzero(numpy.isin(nucleotide_sequence,['a','t','c','g','x','z','o']))==len(nucleotide_sequence), name+": sequence can only contain a,t,c,g,z,x"

        self.nucleotide_sequence=nucleotide_sequence
        self.index=index
        self.nucleotide_number=nucleotide_number
        self.is_cds=is_cds
        self.is_promoter=is_promoter
        self.is_indel=is_indel
        self.indel_length=indel_length

        #Make appropriate changes to the arrays to encorporate the frame shift
        for shift in ribosomal_shifts:
            self.__duplicate(shift)
        if self.reverse_complement:
            self.nucleotide_sequence=self._complement(self.nucleotide_sequence[::-1])
            self.index=self.index[::-1]
            self.nucleotide_number=self.nucleotide_number[::-1]
            self.is_cds=self.is_cds[::-1]
            self.is_promoter=self.is_promoter[::-1]
            self.is_indel=self.is_indel[::-1]
            self.indel_length=self.indel_length[::-1]
            if self.codes_protein:
                self.triplet_number=numpy.floor_divide(self.nucleotide_number[self.is_cds]+2,3)

        else:
            if self.codes_protein:
                self.triplet_number=numpy.floor_divide(self.nucleotide_number[self.is_cds]+2,3)

        self.total_number_nucleotides=len(nucleotide_sequence)

        if self.codes_protein:
            self._setup_conversion_dicts()
            self._translate_sequence()
    
    def __duplicate(self, index):
        '''Duplicate all indcides of important arrays to add the ribosomal shift

        Args:
            index (int): Gene index to duplicate in all arrays
        '''
        #Convert to gene array index
        index = int(numpy.where(self.nucleotide_number == index)[0])
        #TODO: This needs modifying to be correct...
        print(index)
        print(self.nucleotide_number.tolist())

        #Update the nucelotide_numebers so they include the duplicate
        #Check for promoters before the codons
        first_half = [self.nucleotide_number[i] for i in range(index+1)]
        second_half = [self.nucleotide_number[i] + 1 if self.nucleotide_number[i] > 0 else self.nucleotide_number[i] for i in range(index-1, len(self.nucleotide_number))]
        # first_half = [self.nucleotide_number[i] for i in range(len([i for i in self.is_promoter if i])-1 , index)]
        # second_half = [self.nucleotide_number[i]+1 if self.nucleotide_number[i] > 0 else self.nucleotide_number[i] for i in range(index, len(self.nucleotide_number))]
        self.nucleotide_number = numpy.array(first_half + [index] + second_half)
        #Update all 
        self.nucleotide_sequence = self.__duplicate_index(index, self.nucleotide_sequence)
        self.index = self.__duplicate_index(index, self.index)
        self.is_cds = self.__duplicate_index(index, self.is_cds)
        self.is_promoter = self.__duplicate_index(index, self.is_promoter)
        self.is_indel = self.__duplicate_index(index, self.is_indel)
        self.indel_length = self.__duplicate_index(index, self.indel_length)

        

    def __duplicate_index(self, index, array):
        '''Duplicates an element at a given index and returns the new array

        Args:
            index (int): Index of the array to duplicate
            array (numpy.array): Array of items
        Returns:
            numpy.array: Array with duplicated item
        '''
        first_half = [array[i] for i in range(index)]
        second_half = [array[i] for i in range(index, len(array))]
        return numpy.array(first_half + [array[index]] + second_half)

    def __eq__(self, other):
        '''
        Overloading the equality operator to provide a method for determining if two genes
            are the same
        Args:
            other (gumpy.Gene) : The other gene object to compare against
        Returns:
            bool : Boolean showing equality
        '''
        #Default to true
        check = True
        #Check all fields
        check = check and self.name == other.name
        check = check and numpy.all(self.nucleotide_sequence == other.nucleotide_sequence)
        check = check and numpy.all(self.index == other.index)
        check = check and self.__list_eq(self.nucleotide_number, other.nucleotide_number)
        check = check and numpy.all(self.is_cds == other.is_cds)
        check = check and numpy.all(self.is_promoter == other.is_promoter)
        check = check and numpy.all(self.is_indel == other.is_indel)
        check = check and numpy.all(self.indel_length == other.indel_length)
        check = check and numpy.all(self.reverse_complement == other.reverse_complement)
        check = check and numpy.all(self.codes_protein == other.codes_protein)
        check = check and numpy.all(self.feature_type == other.feature_type)
        if self.codes_protein:
            check = check and numpy.all(self.amino_acid_sequence == other.amino_acid_sequence)
            check = check and numpy.all(self.codons == other.codons)
        # if not check:
        #     print(self.name)
        #     print("name", self.name == other.name)
        #     print("NS", numpy.all(self.nucleotide_sequence == other.nucleotide_sequence))
        #     print(self.nucleotide_sequence, self.nucleotide_sequence.shape)
        #     print(other.nucleotide_sequence, other.nucleotide_sequence.shape)
        #     print("index", numpy.all(self.index == other.index))
        #     print("NN", numpy.all(self.nucleotide_number == other.nucleotide_number))
        #     print("is_cds", numpy.all(self.is_cds == other.is_cds))
        #     print("is_promoter", numpy.all(self.is_promoter == other.is_promoter))
        #     print("is_indel", numpy.all(self.is_indel == other.is_indel))
        #     print("indel_len", numpy.all(self.indel_length == other.indel_length))
        #     print("rev_comp", numpy.all(self.reverse_complement == other.reverse_complement))
        #     print("codes_protein", numpy.all(self.codes_protein == other.codes_protein))
        #     print("feature_type", numpy.all(self.feature_type == other.feature_type))
        #     if self.codes_protein:
        #         print("AAS", numpy.all(self.amino_acid_sequence == other.amino_acid_sequence))
        #         print("codons", numpy.all(self.codons == other.codons))
        #     print()
        return check

    @staticmethod
    def _complement(nucleotides_array):
        """
        Simple private method for returning the complement of an array of bases.

        Note that takes account of HET and NULL calls via z and x, respectively
        """

        complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z', 'o':'o', 'n':'n', 'r':'y', 'y':'r', 's':'w', 'w':'s'}

        complement=[complementary_bases[i] for i in nucleotides_array]

        return(numpy.array(complement))

    def _translate_sequence(self):
        #TODO: Optimise this - may require some extra knowledge on why it is done like this...

        # this will ensure that only amino acids with all three bases present
        unique,counts=numpy.unique(self.triplet_number,return_counts=True)
        self.amino_acid_number=unique[counts==3]
        self.amino_acid_number=self.amino_acid_number.astype(int)

        cds_sequence=self.nucleotide_sequence[self.is_cds]

        number_codons=int(numpy.sum(self.is_cds)/3)

        stacked_codons=cds_sequence.reshape((number_codons,3))

        codons=numpy.char.add(stacked_codons[:,0],stacked_codons[:,1])
        self.codons=numpy.char.add(codons,stacked_codons[:,2])

        # now translate the triplets into amino acids using this new dictionary
        self.amino_acid_sequence=numpy.array([self.codon_to_amino_acid[i] for i in self.codons])

    def _setup_conversion_dicts(self):

        bases = ['t', 'c', 'a', 'g', 'x', 'z', 'o']
        # aminoacids = 'FFLLXZSSSSXZYY!!XZCC!WXZXXXXXXZZZZXZLLLLXZPPPPXZHHQQXZRRRRXZXXXXXXZZZZXZIIIMXZTTTTXZNNKKXZSSRRXZXXXXXXZZZZXZVVVVXZAAAAXZDDEEXZGGGGXZXXXXXXZZZZXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZZZZZXZZZZZXZZZZZXZXXXXXXZZZZXZ'
        # aminoacids = 'FFLLXZOSSSSXZOYY!!XZOCC!WXZOXXXXXXOZZZZXZOOOOOOOOLLLLXZOPPPPXZOHHQQXZORRRRXZOXXXXXXOZZZZXZOOOOOOOOIIIMXZOTTTTXZONNKKXZOSSRRXZOXXXXXXOZZZZXZOOOOOOOOVVVVXZOAAAAXZODDEEXZOGGGGXZOXXXXXXOZZZZXZOOOOOOOOXXXXXXOXXXXXXOXXXXXXOXXXXXXOXXXXXXOXXXXXXOOOOOOOOZZZZXZOZZZZXZOZZZZXZOZZZZXZOXXXXXXOZZZZXZOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
        aminoacids = 'FFLLXZOSSSSXZOYY!!XZOCC!WXZOXXXXXXXZZZZXZOOOOOXOOLLLLXZOPPPPXZOHHQQXZORRRRXZOXXXXXXXZZZZXZOOOOOXOOIIIMXZOTTTTXZONNKKXZOSSRRXZOXXXXXXXZZZZXZOOOOOXOOVVVVXZOAAAAXZODDEEXZOGGGGXZOXXXXXXXZZZZXZOOOOOXOOXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZOZZZZXZOZZZZXZOZZZZXZOXXXXXXXZZZZXZOOOOOXOOOOOOXOOOOOOXOOOOOOXOOOOOOXOOXXXXXXXOOOOXOOOOOOXOO'
        all_codons = numpy.array([a+b+c for a in bases for b in bases for c in bases])
        self.codon_to_amino_acid = dict(zip(all_codons, aminoacids))
        # self.amino_acids_of_codons=numpy.array([self.codon_to_amino_acid[i] for i in all_codons])

    def __repr__(self):

        string_length=5

        output=self.name+" gene\n"
        output+="%i nucleotides" % self.total_number_nucleotides
        if self.codes_protein:
            output+=", codes for protein\n"
        else:
            output+="\n"
        if self.nucleotide_sequence[self.is_promoter].size!=0:
            #Updated to use numpy.printoptions to control the formatting
            #Does return as `['a' 'c' 't' 'g']`` rather than `actg` but fixes formatting issues from old code
            #   such as repeating elements if len(promoter) <= string_length
            with numpy.printoptions(threshold=string_length*2):
                output += str(self.nucleotide_sequence[self.is_promoter]) + "\n"
                output += str(self.nucleotide_number[self.is_promoter]) + "\n"
        else:
            output+="promoter likely in adjacent gene(s)\n"
        if self.codes_protein:
            with numpy.printoptions(threshold=string_length*2):
                output += str(self.amino_acid_sequence) + "\n"
                output += str(self.amino_acid_number)
        else:
            with numpy.printoptions(threshold=string_length*2):
                output += str(self.nucleotide_sequence[self.is_cds]) + "\n"
                output += str(self.nucleotide_number[self.is_cds])

        if output.strip()=="":
            output=None

        return(output)

    def list_mutations_wrt(self,other):
        '''Generate a list of mutations between this Gene object and a provided other.
        The mutations are given in the GARC

        Args:
            other (gumpy.Gene): Other Gene object.
        Returns:
            list: List of strings showing the mutations between the two Genes.
        '''

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.name==other.name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        mutations=[]
        if self.codes_protein:

            # deal first with the coding sequence, which we will view at the amino acid level
            # using codons rather than amino acids will detect synonymous mutations as well
            mask=(self.codons!=other.codons)
            pos=self.amino_acid_number[mask]
            ref=other.amino_acid_sequence[mask]
            alt=self.amino_acid_sequence[mask]

            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)

            mask=(self.nucleotide_sequence!=other.nucleotide_sequence) & self.is_promoter
            pos=self.nucleotide_number[mask]
            ref=other.nucleotide_sequence[mask]
            alt=self.nucleotide_sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)
        else:
            mask=self.nucleotide_sequence!=other.nucleotide_sequence
            pos=self.nucleotide_number[mask]
            ref=other.nucleotide_sequence[mask]
            alt=self.nucleotide_sequence[mask]
            for (r,p,a) in zip(ref,pos,alt):
                mutations.append(r+str(int(p))+a)

        mask=self.is_indel
        pos=list(self.nucleotide_number[mask])
        length=list(self.indel_length[mask])
        for (p,l) in zip(pos,length):
            if not l:
                mutations.append(str(int(p))+"_indel_"+str(l))
            elif l > 0:
                mutations.append(str(p)+"_ins_"+str(l))
            else:
                mutations.append(str(p)+"_del_"+str(l))

        if not mutations:
            mutations=None

        return(mutations)
    
    def __list_eq(self, l1, l2):
        '''
        Function used to test equality of 2 lists easily (used in self.__eq__)
        Args:
            l1 (list) : List 1
            l2 (list) : List 2
        Returns:
            bool : Boolean showing whether the 2 lists are element-wise equal
        '''
        if len(l1) != len(l2):
            return False
        for (x1, x2) in zip(l1, l2):
            if x1 != x2:
                return False
        return True



    # def table_mutations_wrt(self,other):
    #
    #     assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
    #     assert self.name==other.name, "both genes must be identical!"
    #     assert self.codes_protein==other.codes_protein, "both genes must be identical!"
    #
    #     MUTATIONS_dict={}
    #     MUTATIONS_columns=['GENE','MUTATION','REF','ALT','POSITION','AMINO_ACID_NUMBER','GENOME_INDEX','NUCLEOTIDE_NUMBER','IS_SNP','IS_INDEL','IS_NULL','IS_FILTER_PASS','IN_CDS','IN_PROMOTER','ELEMENT_TYPE','MUTATION_TYPE','INDEL_LENGTH','INDEL_1','INDEL_2']
    #     for cols in MUTATIONS_columns:
    #         MUTATIONS_dict[cols]=[]
    #
    #     if self.codes_protein:
    #
    #         # deal first with the coding sequence, which we will view at the amino acid level
    #         # using codons rather than amino acids will detect synonymous mutations as well
    #         mask=(self.codons!=other.codons)
    #         pos=self.amino_acid_numbering[mask]
    #         ref=other.amino_acid_sequence[mask]
    #         alt=self.amino_acid_sequence[mask]
    #         ref_codon=other.codons[mask]
    #         alt_codon=self.codons[mask]
    #
    #         for (r,p,a,rc,ac) in zip(ref,pos,alt,ref_codon,alt_codon):
    #             mut=r+str(int(p))+a
    #             MUTATIONS_dict['GENE'].append(self.name)
    #             MUTATIONS_dict['MUTATION'].append(mut)
    #             MUTATIONS_dict['REF'].append(rc)
    #             MUTATIONS_dict['ALT'].append(ac)
    #             MUTATIONS_dict['POSITION'].append(p)
    #             MUTATIONS_dict['AMINO_ACID_NUMBER'].append(p)
    #             MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(0)
    #             MUTATIONS_dict['GENOME_INDEX'].append(0)
    #             MUTATIONS_dict['IS_SNP'].append(True)
    #             MUTATIONS_dict['IS_INDEL'].append(False)
    #             if a=='X':
    #                 MUTATIONS_dict['IS_NULL'].append(True)
    #             else:
    #                 MUTATIONS_dict['IS_NULL'].append(False)
    #             if a=='O':
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(False)
    #             else:
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(True)
    #             MUTATIONS_dict['IN_CDS'].append(True)
    #             MUTATIONS_dict['IN_PROMOTER'].append(False)
    #             MUTATIONS_dict['INDEL_LENGTH'].append(0)
    #             MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
    #             MUTATIONS_dict['MUTATION_TYPE'].append('AAM')
    #             MUTATIONS_dict['INDEL_1'].append('')
    #             MUTATIONS_dict['INDEL_2'].append('')
    #
    #         mask=(self.sequence!=other.sequence) & self.is_promoter
    #         pos=self.numbering[mask]
    #         ref=other.sequence[mask]
    #         alt=self.sequence[mask]
    #         idx=self.index[mask]
    #         for (r,p,a,i) in zip(ref,pos,alt,idx):
    #             mut=r+str(int(p))+a
    #             MUTATIONS_dict['GENE'].append(self.name)
    #             MUTATIONS_dict['MUTATION'].append(mut)
    #             MUTATIONS_dict['REF'].append(r)
    #             MUTATIONS_dict['ALT'].append(a)
    #             MUTATIONS_dict['POSITION'].append(p)
    #             MUTATIONS_dict['AMINO_ACID_NUMBER'].append(0)
    #             MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
    #             MUTATIONS_dict['GENOME_INDEX'].append(i)
    #             MUTATIONS_dict['IS_SNP'].append(True)
    #             MUTATIONS_dict['IS_INDEL'].append(False)
    #             if a=='x':
    #                 MUTATIONS_dict['IS_NULL'].append(True)
    #             else:
    #                 MUTATIONS_dict['IS_NULL'].append(False)
    #             if a=='o':
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(False)
    #             else:
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(True)
    #             MUTATIONS_dict['IN_CDS'].append(False)
    #             MUTATIONS_dict['IN_PROMOTER'].append(True)
    #             MUTATIONS_dict['INDEL_LENGTH'].append(0)
    #             MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
    #             MUTATIONS_dict['MUTATION_TYPE'].append('SNP')
    #             MUTATIONS_dict['INDEL_1'].append('')
    #             MUTATIONS_dict['INDEL_2'].append('')
    #
    #     else:
    #         mask=self.sequence!=other.sequence
    #         pos=self.positions[mask]
    #         ref=other.sequence[mask]
    #         alt=self.sequence[mask]
    #         idx=self.index[mask]
    #         # filter_fail=self.is_filter_fail[mask]
    #         for (r,p,a,i) in zip(ref,pos,alt,idx):
    #             if p<0:
    #                 is_promoter=True
    #                 is_cds=False
    #             else:
    #                 is_promoter=False
    #                 is_cds=True
    #             mut=r+str(int(p))+a
    #             MUTATIONS_dict['GENE'].append(self.name)
    #             MUTATIONS_dict['MUTATION'].append(mut)
    #             MUTATIONS_dict['REF'].append(r)
    #             MUTATIONS_dict['ALT'].append(a)
    #             MUTATIONS_dict['POSITION'].append(p)
    #             MUTATIONS_dict['AMINO_ACID_NUMBER'].append(0)
    #             MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
    #             MUTATIONS_dict['GENOME_INDEX'].append(i)
    #             MUTATIONS_dict['IS_SNP'].append(True)
    #             MUTATIONS_dict['IS_INDEL'].append(False)
    #             if a=='x':
    #                 MUTATIONS_dict['IS_NULL'].append(True)
    #             else:
    #                 MUTATIONS_dict['IS_NULL'].append(False)
    #             if a=='o':
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(False)
    #             else:
    #                 MUTATIONS_dict['IS_FILTER_PASS'].append(True)
    #             MUTATIONS_dict['IN_CDS'].append(is_cds)
    #             MUTATIONS_dict['IN_PROMOTER'].append(is_promoter)
    #             MUTATIONS_dict['INDEL_LENGTH'].append(0)
    #             MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
    #             MUTATIONS_dict['MUTATION_TYPE'].append('SNP')
    #             MUTATIONS_dict['INDEL_1'].append('')
    #             MUTATIONS_dict['INDEL_2'].append('')
    #
    #
    #     mask=self.is_indel
    #     pos=list(self.positions[mask])
    #     num=list(self.numbering[mask])
    #     length=list(self.indel_length[mask])
    #     idx=self.index[mask]
    #     for (p,l,n,i) in zip(pos,length,num,idx):
    #         if p<0:
    #             is_promoter=True
    #             is_cds=False
    #             n=None
    #         else:
    #             is_promoter=False
    #             is_cds=True
    #         mut0=str(int(p))+"_indel"
    #         if l>0:
    #             mut1=str(int(p))+"_ins"
    #             mut2=mut1+"_"+str(l)
    #         else:
    #             mut1=str(int(p))+"_del"
    #             mut2=mut1+"_"+str(-1*l)
    #         MUTATIONS_dict['GENE'].append(self.name)
    #         MUTATIONS_dict['MUTATION'].append(mut0)
    #         MUTATIONS_dict['REF'].append('')
    #         MUTATIONS_dict['ALT'].append('')
    #         MUTATIONS_dict['POSITION'].append(p)
    #         MUTATIONS_dict['AMINO_ACID_NUMBER'].append(n)
    #         MUTATIONS_dict['NUCLEOTIDE_NUMBER'].append(p)
    #         MUTATIONS_dict['GENOME_INDEX'].append(i)
    #         MUTATIONS_dict['IS_SNP'].append(False)
    #         MUTATIONS_dict['IS_INDEL'].append(True)
    #         MUTATIONS_dict['IS_NULL'].append(False)
    #         MUTATIONS_dict['IS_FILTER_PASS'].append(True)
    #         MUTATIONS_dict['IN_CDS'].append(is_cds)
    #         MUTATIONS_dict['IN_PROMOTER'].append(is_promoter)
    #         MUTATIONS_dict['INDEL_LENGTH'].append(l)
    #         MUTATIONS_dict['ELEMENT_TYPE'].append(self.gene_type)
    #         MUTATIONS_dict['MUTATION_TYPE'].append('INDEL')
    #         MUTATIONS_dict['INDEL_1'].append(mut1)
    #         MUTATIONS_dict['INDEL_2'].append(mut2)
    #
    #     if len(MUTATIONS_dict['POSITION'])>0:
    #
    #         MUTATIONS_table=pandas.DataFrame(data=MUTATIONS_dict)
    #
    #         MUTATIONS_table=MUTATIONS_table[MUTATIONS_columns]
    #
    #         MUTATIONS_table=MUTATIONS_table.astype({'GENE':'category',\
    #                                                 'MUTATION':'str',\
    #                                                 'REF':'str',\
    #                                                 'ALT':'str',\
    #                                                 'POSITION':'float',\
    #                                                 'AMINO_ACID_NUMBER':'float',\
    #                                                 'NUCLEOTIDE_NUMBER':'float',\
    #                                                 'GENOME_INDEX':'float',\
    #                                                 'IS_SNP':'bool',\
    #                                                 'IS_INDEL':'bool',\
    #                                                 'IS_NULL':'bool',\
    #                                                 'IS_FILTER_PASS':'bool',\
    #                                                 'IN_CDS':'bool',\
    #                                                 'IN_PROMOTER':'bool',\
    #                                                 'INDEL_LENGTH':'float',\
    #                                                 'ELEMENT_TYPE':'category',\
    #                                                 'MUTATION_TYPE':'category',\
    #                                                 'INDEL_1':'str',\
    #                                                 'INDEL_2':'str'})
    #
    #         MUTATIONS_table=MUTATIONS_table.replace({   'POSITION':0,\
    #                                                     'NUCLEOTIDE_NUMBER':0,\
    #                                                     'AMINO_ACID_NUMBER':0,\
    #                                                     'GENOME_INDEX':0,\
    #                                                     'INDEL_LENGTH':0  }, numpy.nan)
    #
    #         return(MUTATIONS_table)
    #
    #     else:
    #
    #         return(None)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genes.
        Differences are given in the form of an array of gene indices where the two Genes differ.
        
        Args:
            other (gumpy.Gene): Other gene object
        Raises:
            AssertationErrors: Raises errors if the Genes are not the same gene (same length, name and codes_protein are required)
        Returns:
            numpy.array: Numpy array of gene indices where the Genes are different.
        """

        assert self.total_number_nucleotides==other.total_number_nucleotides, "genes must have the same length!"
        assert self.name==other.name, "both genes must be identical!"
        assert self.codes_protein==other.codes_protein, "both genes must be identical!"

        positions=[]

        mask=self.nucleotide_sequence!=other.nucleotide_sequence
        positions=list(self.nucleotide_number[mask])

        mask=self.is_indel
        positions+=list(self.nucleotide_number[mask])


        if not positions:
            return None

        return(numpy.array(positions))


    def valid_variant(self, variant):
        '''Determines if a given variant is valid

        Args:
            variant (str): String of a mutation in GARC

        Raises:
            TypeError: TypeError is raised if the variant is malformed
            AssertationError: AssertationError is raised if the vairant is not valid

        Returns:
            bool: Returns True when the variant is valid. Never returns False.
        '''
        #TODO: Fix or remove this?? self.positions does not exist... Is it referring to self.index?
        #Surely it would make more sense to actually return False at some point instead of raising an AssertaionError??

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
        assert before==self.nucleotide_sequence[mask][0], "base in genome is "+self.nucleotide_sequence[mask][0]+" but specified base is "+before

        # check that the base to be mutated to is valid (z=het, ?=any other base according to the grammar)
        assert after in ['c','t','g','a','?','z','x','o'], after+" is not a valid nucleotide!"

        return True
