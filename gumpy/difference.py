import numpy
import warnings
from collections import defaultdict
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
            self.indels = self._indels_full
            if isinstance(self, GeneDifference) and self._codes_protein:
                #Gene specific attributes
                self.codons = self._codons_full
                self.amino_acids = self._amino_acids_full
        if method == "diff":
            #Convert the full arrays into diff arrays
            self.nucleotides = self.__full_to_diff(self._nucleotides_full)
            self.indels = self.__full_to_diff(self._indels_full)
            if isinstance(self, GeneDifference) and self._codes_protein:
                #Gene specific attributes
                self.codons = self.__full_to_diff(self._codons_full)
                self.amino_acids = self.__full_to_diff(self._amino_acids_full)
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

        #Calculate SNPs
        self.snp_distance = self.__snp_distance()
        '''
        Where applicable, the `full` diference arrays are stored as these can be easily converted into `diff` arrays but not the other way around.
        '''
        #Calculate differences
        self.indices = self.__indices()
        self._nucleotides_full = self.__nucleotides()

        #These are only valid when a VCF has been applied to at least 1 of the genomes
        self.indel_indices = self.__indel_indices()
        self._indels_full = self.__indels()

        #Checking for the same genes, give a warning is the genes are different
        if self.genome1.genes.keys() != self.genome2.genes.keys():
            self.__raise_mutations_warning(self.genome1, self.genome2)

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
    Reports the values of indel calls which differ in length from any indels within the genome.
    Whenever an index is given within an array/dict, it refers to the genome index (i.e 1 indexed)
    '''
    def __init__(self, vcf, genome):
        '''VCF difference object constructor.

        Args:
            vcf (gumpy.VariantFile): VariantFile object
            genome (gumpy.Genome): Genome object
        '''
        self.genome = genome
        self.vcf = vcf

        self.__get_variants()
        self.snps = self.__snps()
        self.snp_distance = len(self.snps)

        self.indels = self.__indels()

        self.genes= genome.stacked_gene_name[numpy.isin(genome.stacked_nucleotide_index,(self.indices))]

    def __get_variants(self):
        '''Pull the variants out of the VariantFile object. Builds arrays
            of the variant calls, and their respective genome indices, as well as
            masks to show whether there is a snp, het, null or indel call at the corresponding genome index:
            i.e is_snp[genome.nucleotide_number == indices[i]] gives a bool to determine if a genome has a SNP call at this position
        '''
        calls = []
        indices = []
        is_snp = []
        is_het = []
        is_null = []
        is_indel = []
        variants = defaultdict(list)

        for index in sorted(list(self.vcf.variants.keys())):
            indices.append(index)
            call = self.vcf.variants[index]['call']
            #Update the masks with the appropriate types
            if self.vcf.variants[index]["type"] == 'indel':
                #Convert to ins_x or del_x rather than tuple
                call = call[0]+"_"+str(call[1])
                is_indel.append(True)
                is_snp.append(False)
                is_het.append(False)
                is_null.append(False)
            elif self.vcf.variants[index]["type"] == "snp":
                is_indel.append(False)
                is_snp.append(True)
                is_het.append(False)
                is_null.append(False)
            elif self.vcf.variants[index]['type'] == 'het':
                is_indel.append(False)
                is_snp.append(False)
                is_het.append(True)
                is_null.append(False)
            elif self.vcf.variants[index]['type'] == 'null':
                is_indel.append(False)
                is_snp.append(False)
                is_het.append(False)
                is_null.append(True)
            calls.append(call)
            for key in self.vcf.variants[index]['original_vcf_row']:
                variants[key].append(self.vcf.variants[index]['original_vcf_row'][key])
        #Convert to numpy arrays for neat indexing
        self.calls = numpy.array(calls)
        self.indices = numpy.array(indices)
        self.is_indel = numpy.array(is_indel)
        self.is_snp = numpy.array(is_snp)
        self.is_het = numpy.array(is_het)
        self.is_null = numpy.array(is_null)
        self.variants = dict()
        for key in variants:
            self.variants[key] = numpy.array(variants[key], dtype=object)

    def __snps(self):
        '''Find the SNPs positions caused by this VCF. Het and null calls are included in SNPs.

        Returns:
            dict: Dictionary mapping genome_index->snp_call
        '''
        #Utilise the masks to determine SNP positions
        #Include HET and NULL calls in SNPs too
        mask = numpy.logical_or(
                            numpy.logical_or(self.is_snp, self.is_het),
                            self.is_null)
        #Get dict mapping genome_index->snp_call
        _snps = dict(zip(self.indices[mask],self.calls[mask]))

        #Convert to GARC mutation nomenclature of ref>call
        snps = {}
        for index in _snps.keys():
            snps[index] = self.genome.nucleotide_sequence[index-1]+">"+_snps[index]
        return snps

    def __indels(self):
        '''Find the difference in the indels and the positions which it varies at.

        Returns:
            dict: Dictionary mapping genome_index->indel
        '''
        indels = dict(zip(self.indices[self.is_indel], self.calls[self.is_indel]))
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
        if gene1.codes_protein != gene2.codes_protein:
            warnings.warn(f"The two genes given do not have the same protein coding for {gene1.name}: Gene1 = {gene1.codes_protein}, Gene2 = {gene2.codes_protein}", UserWarning)
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
