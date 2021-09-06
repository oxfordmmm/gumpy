import numpy, warnings
from abc import ABC #Python library for abstract classes

'''
Classes for difference objects
'''

class Difference(ABC):
    '''Abstract class used to provide the ability to switch between views. Inherited by GenomeDifference and GeneDifference.
    This should not be instantiated.
    '''
    def __init__(self):
        #If this class is instantiated, crash
        assert False, "This class should not be instantiated!"
    def update_view(self, method):
        '''Update the viewing method. Can either be `diff` or `full`:
        `diff`: Where applicable, variables return arrays of values from object1 where they are not equal to object 2
        `full`: Where applicable, variables return arrays of tuples showing (object1_val, object2_val) where values are not equal between objects.

        Args:
            method (str): Name of the viewing method. Must be within `['diff', 'full']`
        '''
        assert type(method) == str, "Invalid method "+str(method)+" of type "+str(type(method))
        assert method in ['diff', 'full'], "Invalid method: "+method

        if method == "full":
            self.nucleotides = self._nucleotides_full
            if hasattr(self,'_codes_protein') and self._codes_protein:
                self.codons = self._codons_full
                self.amino_acids = self._amino_acids_full
            self.indels = self._indels_full
            if hasattr(self, "_het_calls_full"):
                self.het_calls = self._het_calls_full
        if method == "diff":
            #Convert the full arrays into diff arrays
            self.nucleotides = self.__full_to_diff(self._nucleotides_full)
            if hasattr(self,'_codes_protein') and self._codes_protein:
                self.codons = self.__full_to_diff(self._codons_full)
                self.amino_acids = self.__full_to_diff(self._amino_acids_full)
            self.indels = self.__full_to_diff(self._indels_full)
            if hasattr(self, "_het_calls_full"):
                self.het_calls = self.__full_to_diff(self._het_calls_full)
        self._view_method = method

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
        if type(arr1) not in [list, tuple, type(numpy.array([]))]:
            return check or (arr1 != arr2)
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
            numpy.array: Array of values from object1
        '''
        if len(array) > 0:
            array = numpy.array([item1 for (item1, item2) in array if self.__check_any(item1, item2, False)], dtype=object)
        else:
            return numpy.array([])
        #Check for all None values
        if self.__check_none(array, True):
            return None
        else:
            return array

class GenomeDifference(Difference):
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
        '''Constructor for the GenomeDifference object. Called implictly when `genome1.difference(genome2)` is performed.
        Args:
            genome1 (gumpy.Genome): A Genome object to compare against
            genome2 (gumpy.Genome): The other Genome object
        '''
        self.genome1 = genome1
        self.genome2 = genome2
        self._view_method = "diff"
        self._codes_protein=False

        #Calculate SNPs
        self.snp_distance = self.__snp_distance()
        '''
        Where applicable, the `full` diference arrays are stored as these can be easily converted into `diff` arrays but not the other way around.
        '''
        #Calculate differences
        self.indices = self.__indices()
        self._nucleotides_full = self.__nucleotides()

        # self._codons_full = self.__codons()
        # self._amino_acids_full = self.__amino_acids()

        #These are only valid when a VCF has been applied to at least 1 of the genomes
        self.indel_indices = self.__indel_indices()
        self._indels_full = self.__indels()

        self.het_indices = self.__het_indices()
        self._het_calls_full = self.__het_calls()

        #Find the mutations if one of the genomes is a reference genome
        #This does not require any re-formatting to meet the `diff` view
        self.mutations = self.__mutations()

        self.update_view(self._view_method)

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
        return self.genome1.nucleotide_index[mask]

    def __nucleotides(self):
        '''Calculate the difference in nucleotides
        Returns:
            numpy.array: Numpy array of tuples of (genome1_nucleotide, genome2_nucleotide)
        '''
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence
        return numpy.array(list(zip(self.genome1.nucleotide_sequence[mask], self.genome2.nucleotide_sequence[mask])))

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
            return numpy.array([(None, indel) for indel in self.genome2.indels.values()], dtype=object)
        elif self.genome2.indels is None:
            return numpy.array([(indel, None) for indel in self.genome1.indels.values()], dtype=object)
        else:
            return numpy.array([(self.genome1.indels.get(index), self.genome2.indels.get(index)) for index in self.indel_indices])

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
            #There are genes in the mutant which are not in the reference
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
        '''Pad lists of mutations to be the same length so zip() doesn't lose results

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
        if self._view_method == "full":
            return numpy.array([
                (m1, m2) for (m1, m2) in
                    list(zip(
                        *self.__pad_mutations(
                            sorted(self.genome1_mutations), sorted(self.genome2_mutations))
                        ))
                if m1 != m2])
        if self._view_method == "diff":
            return sorted(list(set(self.genome1_mutations).difference(set(self.genome2_mutations))))

    def gene_differences(self):
        '''Get the GeneDifference objects for each gene in the genomes.
        Must be explicitly called by the user

        Returns:
            numpy.array: Array of GeneDifference objects
        '''
        #Build the genome with the VCF applied
        #Checking for the same genes
        if self.genome1.genes.keys() != self.genome2.genes.keys():
            #Get only the genes which are the same but give a warning
            genes = set(self.genome1.genes.keys()).intersection(set(self.genome2.genes.keys()))
            self.__raise_mutations_warning(self.genome1, self.genome2)
        else:
            genes = self.genome1.genes.keys()
        return numpy.array([GeneDifference(self.genome1.genes[gene], self.genome2.genes[gene]) for gene in genes])


class VCFDifference(object):
    '''Object used to find the difference a VCF file makes to a given Genome.
    This includes differences within codons/amino acids, coverages and mutations.
    Reports the values of indel calls which differ in length from any indels within the genome
    '''
    def __init__(self, vcf, genome):
        self.genome = genome
        self.vcf = vcf

        self.indices = self.__indices()
        self.snps = self.__snps()
        self.snp_distance = len(self.indices)

        # self.coverages = self.__coverages()
        # self.het_calls = self.__het_calls()
        self.indels = self.__indels()

        self.genes= genome.stacked_gene_name[numpy.isin(genome.stacked_nucleotide_index,(self.indices))]

    def __indices(self):
        '''Find the SNP positions caused by this VCF

        Returns:
            numpy.array: Array of SNP genome indices
        '''
        indices=[]
        for idx in self.vcf.variants:
            if 'indel' not in self.vcf.variants[idx]:
                call=self.vcf.variants[idx]['call']
                if self.genome.nucleotide_sequence[idx-1] != call:
                    indices.append(self.genome.nucleotide_index[idx- 1])
        return numpy.array(indices)

    def __snps(self):
        '''Find the SNPs positions caused by this VCF

        Returns:
            numpy.array: Array of SNP genome indices
        '''
        snps = {}

        for idx in self.vcf.variants:
            # snp, null, het all ok
            if self.vcf.variants[idx]['type']!='indel':
                call=self.vcf.variants[idx]['call']
                if self.genome.nucleotide_sequence[idx-1] != call:
                    snps[idx]=self.genome.nucleotide_sequence[idx- 1]+">"+call
        return numpy.array(snps)

    # def __coverages(self):
    #     '''Finds the coverages of each call at each position
    #
    #     Returns:
    #         dict: Dictionary mapping the genome_index->[(cov, call)]
    #     '''
    #     coverages = {}
    #     for record in self.vcf.records:
    #         if record.values["COV"] == (0, 0) or record.values["COV"] is None:
    #             continue
    #         if record.alts is None:
    #             #Checking for null values
    #             record.alts = ('x', )
    #         coverages[record.pos] = list(zip(record.values["COV"][1::], record.alts))
    #     return coverages
    #
    # def __het_calls(self):
    #     '''Find the possible values for het calls defined in the VCF
    #
    #     Returns:
    #         dict: Dictionary mapping genome_index->[call1, call2..]
    #     '''
    #     het_calls = {}
    #     for record in self.vcf.records:
    #         if len(record.alts) > 1:
    #             #There is a het call
    #             if record.alts is None:
    #                 record.alts = ('x', )
    #             het_calls[record.pos] = record.alts
    #     return het_calls

    def __indels(self):
        '''Find the difference in the indels and the positons which it varies at.

        Returns:
            dict: Dictionary mapping array_index->array(indel)
        '''
        indels = {}
        for index in self.vcf.variants.keys():
            if self.vcf.variants[index]['type']=='indel':
                call=self.vcf.variants[index]['call']
                # Indel call so check for differences
                if self.genome.is_indel[index-1] == True and self.genome.indel_length[index-1] == call[1]:
                    #Check genome.indels to see if they are the same indel

                    if self.genome.indels is not None and numpy.any(self.genome.indels[index-1] != call[0]+"_"+str(call[1])):
                        #There is a same length, but different value indel
                        indels[index] = call[0]+"_"+str(call[1])
                    else:
                        #No differences in indels so ignore it
                        continue
                else:
                    indels[index] = call[0]+"_"+str(call[1])
        return indels


    def gene_differences(self):
        '''Get the GeneDifference objects for each gene in the genome comparing existing genes with genes after VCF.
        Must be explicitly called by the user as applying a VCF is computationally expensive

        Returns:
            numpy.array: Array of GeneDifference objects
        '''
        #Build the genome with the VCF applied
        genome_ = self.genome.apply_variant_file(self.vcf)
        return numpy.array([
            GeneDifference(
                self.genome.genes[gene], genome_.genes[gene]
            )
            for gene in self.genome.genes.keys()
        ])


class GeneDifference(Difference):
    '''Object to store the differences within genes. The same view system is inherited from the Difference class.
    Set to `full` for arrays of tuple values where applicable.
    Set to `diff` for arrays of values from gene1 where the values vary. Default.
    '''
    def __init__(self, gene1, gene2):
        '''Constructor. Takes in two gene objects and calculates the difference in a few areas such as nucleotides, codons, and amino acids.

        Args:
            gene1 (gumpy.Gene): Gene object 1
            gene2 (gumpy.Gene): Gene object 2
        '''
        if gene1.total_number_nucleotides != gene2.total_number_nucleotides:
            #The lengths of the genes are different so comparing them is meaningless
            warnings.warn("The two genes ("+gene1.name+" and "+gene2.name+") are different lengths, so comparision failed...", UserWarning)
            return None
        if gene1.name != gene2.name:
            warnings.warn("The two genes given have different names ("+gene1.name+", "+gene2.name+") but the same length, continuing...", UserWarning)
        self.gene1 = gene1
        self.gene2 = gene2
        self._codes_protein=gene1.codes_protein
        self._view_method = "diff"

        self._nucleotides_full = self.__nucleotides()
        self.mutations = self.__mutations()

        self.indel_indices = self.__indel_indices()
        self._indels_full = self.__indels()
        self._codons_full = self.__codons()
        self._amino_acids_full = self.__amino_acids()

        self.update_view("diff") #Use inherited method to set the view

    def __nucleotides(self):
        '''Find the differences in nucleotides

        Returns:
            numpy.array: Array of tuples (gene1_nucleotide, gene2_nucleotide)
        '''
        return numpy.array(
            [(n1, n2)
                for (n1, n2)
                in zip(self.gene1.nucleotide_sequence, self.gene2.nucleotide_sequence)
                if n1 != n2
            ])
    def __mutations(self):
        '''Generate the list of mutations between the two genes

        Returns:
            numpy.array: Array of mutations in GARC
        '''
        mutations = []
        gene_mutation = self.gene2.list_mutations_wrt(self.gene1)
        if gene_mutation is not None:
            for mutation in gene_mutation:
                mutations.append(self.gene1.name+"@"+mutation)
        return numpy.array(sorted(mutations))
    def __indel_indices(self):
        '''Find the positions at which the indels differ between the two genes

        Returns:
            numpy.array: Array of gene indices where the two genes have different indels
        '''
        mask = self.gene1.indel_length != self.gene2.indel_length
        return self.gene1.nucleotide_number[mask]
    def __indels(self):
        '''Find the lengths of the indels at each position where the two genes' indels differ

        Returns:
            numpy.array: Array of lengths of indels in the form [(gene1_indel, gene2_indel)]
        '''
        mask = self.gene1.indel_length != self.gene2.indel_length
        return numpy.array([
            (i1, i2)
            for (i1, i2) in zip(self.gene1.indel_length[mask], self.gene2.indel_length[mask])
        ])
    def __codons(self):
        '''Find the codon positions which are different within the genes (within codon regions)

        Returns:
            numpy.array: Array of codons which differ of the form [(gene1_codon, gene2_codon)]
        '''
        codons1 = convert_nucleotides_codons(self.gene1.nucleotide_sequence[self.gene1.is_cds])
        codons2 = convert_nucleotides_codons(self.gene2.nucleotide_sequence[self.gene2.is_cds])
        mask = codons1 != codons2
        codons1 = codons1[mask]
        codons2 = codons2[mask]
        return numpy.array(list(zip(codons1, codons2)))


    def __amino_acids(self):
        '''Calculate the difference in amino acid sequences (within codon regions)
        Returns:
            numpy.array: Array of tuples showing (amino_acid_1, amino_acid_2)
        '''
        codon_to_amino_acid = setup_codon_aa_dict()

        aa_diff = []
        for (codon1, codon2) in self._codons_full:
            aa1 = codon_to_amino_acid[codon1]
            aa2 = codon_to_amino_acid[codon2]
            if aa1 != aa2:
                aa_diff.append((aa1, aa2))
        return numpy.array(aa_diff)




'''
Helper functions not specific to a class
'''
def convert_nucleotides_codons(nucleotides):
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
