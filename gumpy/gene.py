'''
Gene object
'''
import re

import numpy

from gumpy import GeneDifference


# FIXME: problems with rrs, mfpB
class Gene(object):

    """Gene object that uses underlying numpy arrays"""

    def __init__(self, *args, **kwargs):
        '''Constructor for the Gene object.

        Args:
            name (str, optional): Name of the gene. Defaults to None
            nucleotide_sequence (numpy.array, optional): Numpy array of the nucleotide sequence. Defaults to None
            nucleotide_index (numpy.array, optional): Numpy array of the gene indices. Defaults to None
            nucleotide_number (numpy.array, optional): Numpy array of the gene numbering. Defaults to None
            is_cds (numpy.array, optional): Numpy array to act as a mask for whether given elements are codons. Defaults to None
            is_promoter (numpy.array, optional): Numpy array to act as a mask for whether given elements are promoters. Defaults to None
            is_indel (numpy.array, optional): Numpy array to act as a mask for whether given elements are indels. Defaults to None
            indel_length (numpy.array, optional): Numpy array denoting the lengths of the indels whenever is_indel==True. Defaults to None
            indel_nucleotides (numpy.array, optional): Numpy array describing the nucleotides inserted or deleted whenever is_indel==True. Defaults to None
            reverse_complement (boolean, optional): Boolean showing whether this gene is reverse complement. Defaults to False
            codes_protein (boolean, optional): Boolean showing whether this gene codes a protein. Defaults to True
            feature_type (str, optional): The name of the type of feature that this gene represents. Defaults to None
            ribosomal_shifts (list(int), optional): Indices of repeated bases due to ribosomal frame shifting. Defaults to []
        '''


        #Set the kwargs
        #UPDATE THIS AS REQUIRED
        allowed_kwargs = ['name', 'nucleotide_sequence', 'nucleotide_index', 'nucleotide_number', 'is_cds', 'is_promoter',
                        'is_indel', 'indel_length', 'indel_nucleotides','codes_protein', 'reverse_complement', 'feature_type',
                        'codon_number', 'total_number_nucleotides', 'codon_to_amino_acid', 'amino_acid_number',
                        'codons', 'amino_acid_sequence', 'ribosomal_shifts']
        #If reloading a Gene, just set the attributes and return
        if "reloading" in kwargs.keys():
            seen = set()
            for key in kwargs.keys():
                if key in allowed_kwargs:
                    setattr(self, key, kwargs[key])
                    seen.add(key)
            for key in set(allowed_kwargs).difference(seen):
                #Default values to None if the kwarg has not been passed
                setattr(self, key, None)

            return

        #Set default values based on kwargs
        name = kwargs.get("name")
        nucleotide_sequence = kwargs.get("nucleotide_sequence")
        nucleotide_index = kwargs.get("nucleotide_index")
        nucleotide_number = kwargs.get("nucleotide_number")
        is_cds = kwargs.get("is_cds")
        is_promoter = kwargs.get("is_promoter")
        is_indel = kwargs.get("is_indel")
        indel_length = kwargs.get("indel_length")
        indel_nucleotides = kwargs.get("indel_nucleotides")
        reverse_complement = kwargs.get("reverse_complement", False)
        codes_protein = kwargs.get("codes_protein", True)
        feature_type = kwargs.get("feature_type")
        ribosomal_shifts = kwargs.get("ribosomal_shifts", [])
        #VCF data which is accessed by GeneDifference. Has no direct references here...

        assert name is not None, "must provide a gene name!"
        assert isinstance(name, str)
        self.name=name

        assert isinstance(feature_type, str)
        self.feature_type=feature_type
        assert codes_protein in [True,False], name+": codes_protein must be True or False!"
        self.codes_protein=codes_protein

        assert reverse_complement in [True,False], name+": reverse_complement must be True or False!"
        self.reverse_complement=reverse_complement

        assert isinstance(nucleotide_sequence,numpy.ndarray), name+": sequence of bases must be a Numpy array!"

        assert isinstance(nucleotide_index,numpy.ndarray) and numpy.issubdtype(nucleotide_index.dtype.type,numpy.integer), name+": genome indices must be a Numpy array of integers!"

        assert isinstance(nucleotide_number,numpy.ndarray), name+": gene numbering must be a Numpy array of integers!"

        # check it is a list, and if it is that it is made up of integers and only integers
        assert isinstance(ribosomal_shifts,list)
        if len(ribosomal_shifts)>0:
            assert all([isinstance(i,int) for i in ribosomal_shifts])

        assert len(nucleotide_index)==len(nucleotide_sequence), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(nucleotide_number), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(is_cds), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(is_promoter), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(is_indel), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(indel_length), 'all inputs arrays must be the same length!'
        assert len(nucleotide_index)==len(indel_nucleotides), 'all inputs arrays must be the same length!'

        nucleotide_sequence=numpy.char.lower(nucleotide_sequence)
        assert numpy.count_nonzero(numpy.isin(nucleotide_sequence,['a','t','c','g','x','z','o']))==len(nucleotide_sequence), name+": sequence can only contain a,t,c,g,z,x"

        self.nucleotide_sequence=nucleotide_sequence
        self.nucleotide_index=nucleotide_index
        self.nucleotide_number=nucleotide_number
        self.is_cds=is_cds
        self.is_promoter=is_promoter
        self.is_indel=is_indel
        self.indel_nucleotides=indel_nucleotides
        self.indel_length=indel_length

        #Make appropriate changes to the arrays to encorporate the frame shift
        for shift in ribosomal_shifts:
            self.__duplicate(shift)
        if self.reverse_complement:
            self.nucleotide_sequence=self._complement(self.nucleotide_sequence[::-1])
            self.nucleotide_index=self.nucleotide_index[::-1]
            self.nucleotide_number=self.nucleotide_number[::-1]
            self.is_cds=self.is_cds[::-1]
            self.is_promoter=self.is_promoter[::-1]
            self.is_indel=self.is_indel[::-1]
            self.indel_length=self.indel_length[::-1]
            if self.codes_protein:
                self.codon_number=numpy.floor_divide(self.nucleotide_number[self.is_cds]+2,3)
                self.gene_position=numpy.concatenate((self.nucleotide_number[self.nucleotide_number<0], self.codon_number))
            else:
                self.gene_position=self.nucleotide_number
        else:
            if self.codes_protein:
                self.codon_number=numpy.floor_divide(self.nucleotide_number[self.is_cds]+2,3)
                self.gene_position=numpy.concatenate((self.nucleotide_number[self.nucleotide_number<0], self.codon_number))
            else:
                self.gene_position=self.nucleotide_number

        self.total_number_nucleotides=len(nucleotide_sequence)

        if self.codes_protein:
            self._setup_conversion_dicts()
            self._translate_sequence()

    def __duplicate(self, index):
        '''Duplicate all indices of important arrays to add the ribosomal shift

        Args:
            index (int): Gene index to duplicate in all arrays
        '''
        #Convert to gene array index
        index = int(numpy.where(self.nucleotide_number == index)[0])+1

        #Update the nucelotide_numbers so they include the duplicate
        #Check for promoters before the codons
        first_half = [self.nucleotide_number[i] for i in range(index)]
        second_half = [
            self.nucleotide_number[i] + 1 if self.nucleotide_number[i] > 0
            else self.nucleotide_number[i]
            for i in range(index, len(self.nucleotide_number))]
        self.nucleotide_number = numpy.array(first_half + [self.nucleotide_number[index]] + second_half)
        #Update all
        self.nucleotide_sequence = self.__duplicate_index(index, self.nucleotide_sequence)
        self.nucleotide_index = self.__duplicate_index(index, self.nucleotide_index)
        self.is_cds = self.__duplicate_index(index, self.is_cds)
        self.is_promoter = self.__duplicate_index(index, self.is_promoter)
        self.is_indel = self.__duplicate_index(index, self.is_indel)
        self.indel_length = self.__duplicate_index(index, self.indel_length)
        self.indel_nucleotides = self.__duplicate_index(index, self.indel_nucleotides)


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
        '''Overloading the equality operator to provide a method for determining if two genes
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
        check = check and numpy.all(self.nucleotide_index == other.nucleotide_index)
        check = check and numpy.all(self.nucleotide_number == other.nucleotide_number)
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

        return check

    @staticmethod
    def _complement(nucleotides_array):
        """Simple private method for returning the complement of an array of bases.

        Note that takes account of HET and NULL calls via z and x, respectively
        """

        complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z', 'o':'o', 'n':'n', 'r':'y', 'y':'r', 's':'w', 'w':'s'}

        complement=[complementary_bases[i] for i in nucleotides_array]

        return(numpy.array(complement))

    def _translate_sequence(self):

        # this will ensure that only amino acids with all three bases present
        unique,counts=numpy.unique(self.codon_number,return_counts=True)
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
        '''Overload the print function to write a summary of the Gene.

        Returns:
            str: String describing the Gene
        '''
        string_length=5

        output=self.name+" gene\n"
        output+="%i nucleotides" % self.total_number_nucleotides
        if self.codes_protein:
            output+=", codes for protein\n"
        else:
            output+="\n"
        if self.nucleotide_sequence[self.is_promoter].size!=0:
            # updated to use numpy.printoptions to control the formatting
            # does return as `['a' 'c' 't' 'g']`` rather than `actg` but fixes formatting issues from old code
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

    def __len__(self):
        '''Return the number of nucleotides in the coding region (i.e. ignoring any assumed promoter)

        Returns:
            int
        '''

        return(len(self.nucleotide_sequence[self.is_cds]))

    def __sub__(self, other):
        '''Generate a GeneDifference object for an in-depth examination of the difference between two genes.

        Args:
            other (gumpy.Gene): Other Gene object

        Returns:
            gumpy.GeneDifference: The GeneDifference object to detail changes at all levels such as indel and amino acid.
        '''

        assert isinstance(other,Gene)

        return GeneDifference(self, other)

    def valid_variant(self, variant):
        '''Determines if a given variant is valid for this gene

        Args:
            variant (str): String of a mutation in GARC
            
        Returns:
            bool: True when variant is valid, False otherwise
        '''
        assert variant is not None, "No Variant given!"
        assert type(variant) == str
        assert len(variant) >= 2, "Variant must be at least 2 characters e.g. A="

        #Match mutation formats using regex

        #Match promoter/non-coding SNP format
        promoter = re.compile(r"""
                ([a-zA-Z_0-9]+@)? #Possibly a leading gene name
                ([acgtzx]) #Reference base
                (-?[0-9]+) #Position
                ([acgtzx]) #Alt base
                """, re.VERBOSE)
        if promoter.fullmatch(variant):
            #The variant is either promoter or non-coding
            #Check that the mutation is valid
            name, ref, pos, alt = promoter.fullmatch(variant).groups()
            valid = True
            #Strip '@' from gene name and compare
            if name is not None and name != "":
                valid = valid and name[:-1] == self.name
            #Check that the pos is in the correct range
            valid = valid and int(pos) in self.nucleotide_number
            #Check that the ref matches the seq at the given pos
            valid = valid and self.nucleotide_sequence[self.nucleotide_number == int(pos)] == ref
            #Mutation should not have same ref and alt
            valid = valid and ref != alt
            return valid
        #Match amino acid SNP
        snp = re.compile(r"""
                    ([a-zA-Z_0-9]+@)? #Possibly leading gene name
                    ([A-Z!]) #Reference amino acid
                    ([0-9]+) #Position
                    ([A-Z!]) #Alt amino acid
                    """, re.VERBOSE)
        if self.codes_protein and snp.fullmatch(variant):
            #The variant is an amino acid SNP
            #Check it is valid
            name, ref, pos, alt = snp.fullmatch(variant).groups()
            valid = True
            #Strip '@' from gene name and compare
            if name is not None and name != "":
                valid = valid and name[:-1] == self.name
            #Check pos is in correct range
            valid = valid and int(pos) in self.amino_acid_number
            #Check ref matches the aa seq at the given pos
            valid = valid and self.amino_acid_sequence[self.amino_acid_number == int(pos)] == ref
            return valid

        #Match amino acid synon-mutation
        synon = re.compile(r"""
                        ([a-zA-Z_0-9]+@)? #Possibly leading gene name
                        ([0-9]+) #Pos
                        = #Synonymous amino acid
                        """, re.VERBOSE)
        if self.codes_protein and synon.fullmatch(variant):
            #Variant is a synonymous mutation
            #Check it is valid
            name, pos = synon.fullmatch(variant).groups()
            valid = True
            if name is not None and name != "":
                valid = valid and name[:-1] == self.name
            valid = valid and int(pos) in self.amino_acid_number
            return valid

        #Match indel
        indel = re.compile(r"""
                    ([a-zA-Z_0-9]+@)? #Possibly leading gene name
                    (-?[0-9]+) #Position
                    _(ins|del|indel)_? #Type
                    ([0-9]+|[acgtzx]+)? #Bases deleted/inserted
                    """, re.VERBOSE)
        if indel.fullmatch(variant):
            #Variant is an indel
            #Check it is valid
            name, pos, type_, bases = indel.fullmatch(variant).groups()
            valid = True
            pos = int(pos)
            #Check the names match
            if name is not None and name != "":
                valid = valid and name[:-1] == self.name
            #Check pos in correct range
            valid = valid and int(pos) in self.nucleotide_number
            if type_ == "indel":
                #If a mutation is given as `indel`, no length/bases should be given
                valid = valid and (bases is None or bases == "")
            if type_ == "del" and bases is not None and bases != "":
                #Mutation was del, so check if the bases given match the ref
                if bases.isnumeric():
                    #A length was given rather than bases so just check that all bases are in the correct range
                    if pos < 0:
                        #Is a promoter so be careful about pos+bases not being 0
                        pos -= 1
                    for i in range(1, int(bases)+1):
                        valid = valid and pos + i in self.nucleotide_number
                else:
                    #Bases were given rather than a length, so check for equality against base seq
                    for index in range(len(bases)):
                        valid = valid and int(pos)+index in self.nucleotide_number
                        valid = valid and self.nucleotide_sequence[self.nucleotide_number == int(pos)+index] == bases[index]
            return valid
        #If none of the mutations have matched, return False
        return False
