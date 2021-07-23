import numpy, warnings

'''
Classes for difference objects
'''

class GenomeDifference(object):
    '''GenomeDifference object is used to return the difference between two genomes.
        The difference can be viewed in one of two ways:
            `diff`: Arrays of the values from genome1 where there are different values in genome2, 
                    or values which exist in genome1 but not genome2 depending on which is more appropriate (closer to classical subtraction).
                    This is the default view.
            `full`: Arrays of tuples consisting of (genome1_val, genome2_val) - more information but less like classical subtraction 
                    and more difficult to wrangle meaningful information from.
        This option can be set by using the `update_view()` method

        Instance variables:
            `snp` (int): SNP distance between the two genomes
            `indices` (numpy.array): Array of genome indices where the two genomes differ in nucleotides
            `nucleotides` (numpy.array): Array of differences in nucleotides. Exact format depends on which view is selected.
            `codons` (numpy.array): Array of differences in codons. Exact format depends on which view is selected.
            `amino_acids` (numpy.array): Array of differences in amino acids. Exact format depends on which view is selected.
            `indel_indices` (numpy.array): Array of indices where the two genomes have indels
            `indels` (numpy.array): Array of differences in inels. Exact format depends on which view is selected.
            `het_indices` (numpy.array): Array of indices where there are het calls in either genome
            `het_calls` (numpy.array): Array of het_calls at `het_indices`. Exact format depends on which view is selected.
            `mutations` (numpy.array): Array of mutations within genes of a mutant genome given one of the genomes is a reference genome.
            `genome1_mutations` (numpy.array): Array of mutations within genes of genome1 compared to a reference. Only set once `find_mutations()` has been called
            `genome2_mutations` (numpy.array): Array of mutations within genes of genome2 compared to a reference. Only set once `find_mutations()` has been called
        Methods:
            `update_view(method: str) -> None`: Used to change the viewing method for instance variables. `method` values are either `diff` or `full`
            `find_mutations(reference: gumpy.Genome) -> numpy.array`: Used to find the mutations within the genes of the geomes when neither are reference genomes. Exact format depends on which view is selected.
    '''
    def __init__(self, genome1, genome2):
        '''Constructor for the GenomeDifference object. Called implictly when `genome1 - genome2` is performed.
        Args:
            genome1 (gumpy.Genome): A Genome object to compare against
            genome2 (gumpy.Genome): The other Genome object
        '''
        self.genome1 = genome1
        self.genome2 = genome2
        self.__view_method = "diff"

        #Calculate SNPs
        self.snp = self.__snp_distance()
        '''
        Where applicable, the `full` diference arrays are stored as these can be easily converted into `diff` arrays but not the other way around.
        '''
        #Calculate differences
        self.indices = self.__indices()
        self.__nucleotides_full = self.__nucleotides()

        self.__codons_full = self.__codons()
        self.__amino_acids_full = self.__amino_acids()

        #These are only valid when a VCF has been applied to at least 1 of the genomes
        self.indel_indices = self.__indel_indices()
        self.__indels_full = self.__indels()

        self.het_indices = self.__het_indices()
        self.__het_calls_full = self.__het_calls()

        #Find the mutations if one of the genomes is a reference genome
        #This does not require any re-formatting to meet the `diff` view
        self.mutations = self.__mutations()

        self.update_view(self.__view_method)
    
    def update_view(self, method):
        '''Update the viewing method. Can either be `diff` or `full`:
        `diff`: Where applicable, variables return arrays of values from genome1 where they are not equal to genome 2
        `full`: Where applicable, variables return arrays of tuples showing (genome1_val, genome2_val) where values are not equal between genomes.

        Args:
            method (str): Name of the viewing method. Must be within `['diff', 'full']`
        '''        
        assert type(method) == str, "Invalid method "+str(method)+" of type "+str(type(method))
        assert method in ['diff', 'full'], "Invalid method: "+method
        if method == "full":
            self.nucleotides = self.__nucleotides_full
            self.codons = self.__codons_full
            self.amino_acids = self.__amino_acids_full
            self.indels = self.__indels_full
            self.het_calls = self.__het_calls_full
        if method == "diff":
            #Convert the full arrays into diff arrays
            self.nucleotides = self.__full_to_diff(self.__nucleotides_full)
            self.codons = self.__full_to_diff(self.__codons_full)
            self.amino_acids = self.__full_to_diff(self.__amino_acids_full)
            self.indels = self.__full_to_diff(self.__indels_full)
            self.het_calls = self.__full_to_diff(self.__het_calls_full)
        self.__view_method = method
    
    def __check_none(self, arr, check):
        '''Recursive function to determine if a given array is all None values.
        Due to the shape of some arrays, implicit equality checking causes warnings

        Args:
            arr (array-like): Array to check
            check (bool): Boolean to propagate the check
        Returns:
            bool: True when all None values in all arrays/sub-arrays
        '''
        if check == False:
            #No point in continuing
            return False
        for elem in arr:
            if type(elem) in [list, tuple, type(numpy.array([]))]:
                check = check and self.__check_none(elem, check)
            else:
                check = check and (elem is None)
        return check
    
    def __check_any(self, arr1, arr2, check):
        '''Recursive function to determine if given arrays are the different at any point.
        Required due to implicit equality of weird shape arrays causing warnings

        Args:
            arr1 (array-like): Array 1
            arr2 (array-like): Array 2
            check (bool): Boolean accumulator

        Returns:
            bool: True when the two lists are different at any point
        '''      
        if type(arr1) != type(arr2):
            return True
        for (e1, e2) in zip(arr1, arr2):
            if type(e1) in [list, tuple, type(numpy.array([]))]:
                check = check or self.__check_any(e1, e2, check)
            else:
                check = check or (e1 != e2)
        return check  
    
    def __full_to_diff(self, array):
        '''Convert an array from a full view to a diff view

        Args:
            array (numpy.array): Array of tuples of values
        Returns:
            numpy.array: Array of values from genome1
        '''        
        array = numpy.array([item1 for (item1, item2) in array if self.__check_any(item1, item2, False)], dtype=object)
        #Check for all None values
        if self.__check_none(array, True):
            return None
        else:
            return array

    def __snp_distance(self):
        '''Calculates the SNP distance between the two genomes

        Returns:
            int: The SNP distance between the two genomes
        '''        
        #Counts differences of `z` and `x` as SNPs - FIXME if required
        return sum([1 for (base1, base2) in zip(self.genome1.nucleotide_sequence, self.genome2.nucleotide_sequence) if base1 != base2])
    
    def __indices(self):
        '''Calculates the indices where the two genomes differ in nucleotide sequence

        Returns:
            numpy.array: Numpy array for the difference in nucleotide sequence
        '''        
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence
        mask += self.genome1.indel_length!=0
        return self.genome1.nucleotide_index[mask]
    
    def __nucleotides(self):
        '''Calculate the difference in nucleotides
        Returns:
            numpy.array: Numpy array of tuples of (genome1_nucleotide, genome2_nucleotide)
        '''
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence
        return numpy.array(list(zip(self.genome1.nucleotide_sequence[mask], self.genome2.nucleotide_sequence[mask])))
    
    def __codons(self):
        '''Convert the nucleotide arrays to arrays of codons. Returns codons which are different
        Returns:
            numpy.array: Numpy array of tuples of (codon1, codon2)
        '''        
        if len(self.__nucleotides_full) == 0:
            #No changes in nucleotides so no difference in codons
            return numpy.array([])
        codons1 = self.__convert_nucleotides_codons(self.genome1.nucleotide_sequence)
        codons2 = self.__convert_nucleotides_codons(self.genome2.nucleotide_sequence)
        mask = codons1 != codons2
        codons1 = codons1[mask]
        codons2 = codons2[mask]
        return numpy.array(list(zip(codons1, codons2)))
    
    def __convert_nucleotides_codons(self, nucleotides):
        '''Helper function to convert an array of nucleotides into an array of codons

        Args:
            nucleotides (numpy.array): Array of nucleotides

        Returns:
            numpy.array: Array of codons
        '''        
        codons = []
        c = ""
        for index in range(1, len(nucleotides)+1):
            c += nucleotides[index-1]
            if index % 3 == 0:
                #There have been 3 bases seen so add the codon
                codons.append(c)
                c = ""
        return numpy.array(codons)
    
    def __amino_acids(self):
        '''Calculate the difference in amino acid sequences
        Returns:
            numpy.array: Array of tuples showing (amino_acid_1, amino_acid_2)
        '''
        codon_to_amino_acid = setup_codon_aa_dict()
        
        #Get the difference in codons (if any)
        if len(self.__codons_full) == 0:
            #No different codons so no different amino acids
            return numpy.array([])
        
        aa_diff = []
        for (codon1, codon2) in self.__codons_full:
            aa1 = codon_to_amino_acid[codon1]
            aa2 = codon_to_amino_acid[codon2]
            if aa1 != aa2:
                aa_diff.append((aa1, aa2))
        return numpy.array(aa_diff)
    
    def __indel_indices(self):
        '''Finds the positions at which there are indels in either genome. Use genome.indels.get(item) for safe retrieval of indels

        Returns:
            numpy.array: Array of array indices where there are indels in either genome
        '''
        if self.genome1.indels is None:
            indices1 = set()
        else:
            indices1 = set(self.genome1.indels.keys())

        if self.genome2.indels is None:
            indices2 = set()
        else:
            indices2 = set(self.genome2.indels.keys())

        return numpy.array(sorted(list(indices1.union(indices2))))
    
    def __indels(self):
        '''Find the indels for both genomes and report where they are different
        Returns:
            numpy.array: Array of tuples showing (indel1, indel2), if an indel is given as None, it did not exist in the genome
        '''        
        #Find the indels for both genomes
        #As Genome.indels defaults to None if a VCF has not been applied, checks are required
        if self.genome1.indels is None and self.genome2.indels is None:
            return numpy.array([])
        elif self.genome1.indels is None:
            return numpy.array([(None, indel[0]) for indel in self.genome2.indels.values()], dtype=object)
        elif self.genome2.indels is None:
            return numpy.array([(indel[0], None) for indel in self.genome1.indels.values()], dtype=object)
        else:
            return numpy.array([(self.genome1.indels.get(index)[0], self.genome2.indels.get(index)[0]) for index in self.indel_indices])
    
    def __het_indices(self):
        '''Find the array indices at which there are het calls

        Returns:
            numpy.array: Array of indices where there are het calls in at least one genome
        '''        
        mask = numpy.logical_or(self.genome1.nucleotide_sequence == 'z', self.genome2.nucleotide_sequence == 'z')
        return numpy.array(range(len(self.genome1)))[mask]
    
    def __het_calls(self):
        '''Find the het calls for each het index for each genome. Will only work if a VCF file has been applied to produce het calls.
        Returns:
            numpy.array: Array of tuples (het_calls1, het_calls2) where het_callsx is an array or None
        '''
        if self.genome1.calls is None and self.genome2.calls is None:
            return numpy.array([])
        elif self.genome1.calls is None:
            return numpy.array([(None, self.genome2.calls[index]) for index in self.het_indices], dtype=object)
        elif self.genome2.calls is None:
            return numpy.array([(self.genome1.calls[index], None) for index in self.het_indices], dtype=object)
        else:
            return numpy.array([(self.genome1.calls.get(index), self.genome2.calls.get(index)) for index in self.het_indices], dtype=object)
    
    def __raise_mutations_warning(self, reference, mutant):
        '''Give a warning to the user that the genes within the two genomes are different.
        Warning displays names of the genes which differ.

        Args:
            reference (gumpy.Genome): Genome object for a reference genome
            mutant (gumpy.Genome): Genome object for a mutated genome
        '''        
        message = "The genomes do not have the same genes. "
        if len(set(reference.genes.keys()).difference(set(mutant.genes.keys()))) > 0:
            #There are genes in the reference which are not in the mutant
            message += "The reference has genes ("
            for gene in set(reference.genes.keys()).difference(set(mutant.genes.keys())):
                message += gene+", "
            message = message[:-2]#Remove trailing comma
            message += ") which are not present in the mutant. "
        if len(set(mutant.genes.keys()).difference(set(reference.genes.keys()))) > 0:
            #There are genes in the reference which are not in the mutant
            message += "The mutant has genes ("
            for gene in set(mutant.genes.keys()).difference(set(reference.genes.keys())):
                message += gene+", "
            message = message[:-2]#Remove trailing comma
            message += ") which are not present in the reference. "
        message += "Continuing only with genes which exist in both genomes."
        warnings.warn(message, UserWarning)
    
    def __mutations(self, reference=None, mutant=None):
        '''Find the mutations within genes. Mutations outside of genes are not considered. FIXME if this is required.

        Returns:
            numpy.array: Array of mutations in GARC
        '''        
        if reference is None and mutant is None:
            #Use XOR to determine if there is 1 reference
            if self.genome1.is_reference ^ self.genome2.is_reference:
                #There is exactly 1 reference genome so mutations can be found
                if self.genome1.is_reference == True:
                    reference = self.genome1
                    mutant = self.genome2
                else:
                    reference = self.genome2
                    mutant = self.genome1
            else:
                return numpy.array([])

        #Checking for the same genes
        if reference.genes.keys() != mutant.genes.keys():
            #Get only the genes which are the same but give a warning
            genes = set(reference.genes.keys()).intersection(set(mutant.genes.keys()))
            self.__raise_mutations_warning(reference, mutant)
        else:
            genes = reference.genes.keys()

        mutations = []
        for gene in genes:
            gene_mutation = reference.genes[gene].list_mutations_wrt(mutant.genes[gene])
            if gene_mutation is not None:
                for mutation in gene_mutation:
                    mutations.append(gene+"@"+mutation)
        return numpy.array(sorted(mutations))
    
    def __pad_mutations(self, arr1, arr2):
        '''Pad lists of mutations to be the same length so zip() doesn't loose results

        Args:
            arr1 (numpy.array): Array1
            arr2 (numpy.array): Array2
        Returns:
            numpy.array, numpy.array: The two arrays padded to be the same length
        '''        
        while len(arr1) > len(arr2):
            arr2 = numpy.append(arr2, None)
        while len(arr2) > len(arr1):
            arr1 = numpy.append(arr1, None)
        return arr1, arr2
    
    def find_mutations(self, reference):
        '''Version of self.__mutations() which takes a reference genome for when neither genome is a reference - finding the difference in mutations.

        Args:
            reference (gumpy.Genome): A reference genome

        Returns:
            (numpy.array): Array of mutations. Structure depends on viewing method set. 
                            `diff` (default): returns an array of mutations present in genome1 but not genome2
                            `full`: returns an array of tuples of (genome1_mutation, genome2_mutation)
        '''
        assert reference.is_reference == True, "Genome passed is not a reference genome!"
        self.genome1_mutations = self.__mutations(reference=reference, mutant=self.genome1)
        self.genome2_mutations = self.__mutations(reference=reference, mutant=self.genome2)
        if self.__view_method == "full":
            return numpy.array(
                list(zip(
                    *self.__pad_mutations(
                        sorted(self.genome1_mutations), sorted(self.genome2_mutations))
                    )))
        if self.__view_method == "diff":
            return sorted(list(set(self.genome1_mutations).difference(set(self.genome2_mutations))))



'''
Helper functions not specific to a class
'''

def setup_codon_aa_dict():
    '''Setup a conversion dictionary to convert codons to amino acids

    Returns:
        dict: Dictionary mapping codons of form 'xyz' to amino acids of form 'X'
    '''    
    #Defined bases
    bases = ['t', 'c', 'a', 'g', 'x', 'z', 'o']
    #Defined amino acids in correct order
    aminoacids = 'FFLLXZOSSSSXZOYY!!XZOCC!WXZOXXXXXXXZZZZXZOOOOOXOOLLLLXZOPPPPXZOHHQQXZORRRRXZOXXXXXXXZZZZXZOOOOOXOOIIIMXZOTTTTXZONNKKXZOSSRRXZOXXXXXXXZZZZXZOOOOOXOOVVVVXZOAAAAXZODDEEXZOGGGGXZOXXXXXXXZZZZXZOOOOOXOOXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZOZZZZXZOZZZZXZOZZZZXZOXXXXXXXZZZZXZOOOOOXOOOOOOXOOOOOOXOOOOOOXOOOOOOXOOXXXXXXXOOOOXOOOOOOXOO'
    all_codons = [a+b+c for a in bases for b in bases for c in bases]
    return dict(zip(all_codons, aminoacids))