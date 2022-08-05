'''
Used to find differences between Genomes and Genes, as well as the impact a VCF has on a genome.

Abstract classes:
    Difference - Abstract class providing the ability to change data views
Classes:
    GenomeDifference
    GeneDifference
    GeneticVariation
Functions:
    convert_nucleotides_codons(numpy.array) -> numpy.array: Converts an array of nucleotides to an array of codons.
    setup_codon_aa_dict() -> dict: Returns a dictionary mapping codon->amino_acid
    collapse_inner_dict(dict) -> dict: Converts a dictionary with an inner dictionary to a single dictionary. Key collisions are not considered
                                        intentionally as this should not be an issue with this use case
'''
import warnings
from abc import ABC  # Python library for abstract classes

import numpy

import gumpy


class Difference(ABC):
    '''
    Abstract class used to provide the ability to switch between views. Inherited by GenomeDifference and GeneDifference.

    This should not be instantiated.
    '''

    def __init__(self):
        '''If this class is instantiated, crash
        '''
        assert False, "This class should not be instantiated!"

    def update_view(self, method, object_type):
        '''Update the viewing method. Can either be `diff` or `full`:
        `diff`: Where applicable, variables return arrays of values from object1 where they are not equal to object 2
        `full`: Where applicable, variables return arrays of tuples showing (object1_val, object2_val) where values are not equal between objects.

        Args:
            method (str): Name of the viewing method. Must be within `['diff', 'full']`
        '''
        assert type(method) == str, "Invalid method "+str(method)+" of type "+str(type(method))
        assert method in ['diff', 'full'], "Invalid method: "+method

        if method == "full":
            # self.indels = self._indels_full
            # if object_type=='gene' and self.codes_protein:
            if isinstance(self, GeneDifference):
                if self.codes_protein:
                    self.amino_acid_sequence = self._amino_acids_full
                else:
                    self.nucleotides=self._nucleotides_full
            else:
                self.nucleotides=self._nucleotides_full

        if method == "diff":
            #Convert the full arrays into diff arrays
            self.nucleotides = self.__full_to_diff(self._nucleotides_full)
            # self.indels = self.__full_to_diff(self._indels_full)
            # if object_type=='gene' and  self.codes_protein:
            if isinstance(self, GeneDifference) and self.codes_protein:
                #Gene specific attributes
                # self.codons = self.__full_to_diff(self._codons_full)
                self.amino_acid_sequence = self.__full_to_diff(self._amino_acids_full, amino_acid=True)
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

    def __full_to_diff(self, array, amino_acid=False):
        '''Convert an array from a full view to a diff view

        Args:
            array (numpy.array): Array of tuples of values
            amino_acid (bool, optional): Optional flag to denote amino acids to ensure synonymous mutations are retained. Defaults to False
        Returns:
            numpy.array: Array of values from object1
        '''
        if len(array) > 0:
            if amino_acid:
                #Amino acids in diff should just be item2 consistently
                #The actual diff should have already been applied...
                array = numpy.array([item2 for (item1, item2) in array], dtype=object)
            else:
                array = numpy.array([item2 for (item1, item2) in array if self.__check_any(item1, item2, False)], dtype=object)
        else:
            #Returning an empty array can cause issues with indels in amino acid sequence not producing same length arrays 
            #as they are handled separately to other arrays which may contain [None]. This breaks the idea of 1:1 relationship
            #between values in arrays. So just return None instead of an array
            return None
        #Check for all None values
        if self.__check_none(array, True):
            return None
        else:
            return array

class GenomeDifference(Difference):
    '''
    GenomeDifference object captures the difference between two genomes.

    The difference can be viewed in one of two ways:
        `diff`: Arrays of the values from genome1 where there are different values in genome2,
                or values which exist in genome1 but not genome2 depending on which is more appropriate (closer to classical subtraction).
                This is the default view.
        `full`: Arrays of tuples consisting of (genome1_val, genome2_val) - more information but less like classical subtraction
                and more difficult to wrangle meaningful information from.
    This option can be set by using the update_view() method

    Instance variables:
        snp_distance (int): SNP distance between the two genomes
        indices (numpy.array): Array of genome indices where the two genomes differ in nucleotides
        nucleotides (numpy.array): Array of differences in nucleotides. Format depends on the current view.
        indel_indices (numpy.array): Array of indices where the two genomes have indels
        indels (numpy.array): Array of differences in inels. Format depends on the current view.
    Functions:
        variants(int) -> dict: Takes a genome index and returns a dictionary mapping field->(genome1_val, genome2_val) for all fields
                                of a vcf file (if applicable)
        gene_differences() -> [GeneDifference]: Returns a list of GeneDifference objects
    Inherited functions:
        update_view(str) -> None: Used to change the viewing method for instance variables. Input values are either `diff` or `full`
    '''
    def __init__(self, genome1, genome2):
        '''
        Constructor for the GenomeDifference object.

        Called implictly when `genome1.difference(genome2)` is performed.
        Args:
            genome1 (gumpy.Genome): A Genome object to compare against
            genome2 (gumpy.Genome): The other Genome object
        '''

        # insist that both must be Genome objects
        assert isinstance(genome1, gumpy.Genome)
        assert isinstance(genome2, gumpy.Genome)

        self.genome1 = genome1
        self.genome2 = genome2
        self._view_method = "diff"

        #Calculate SNPs
        self.snp_distance = self.__snp_distance()
        '''
        Where applicable, the `full` difference arrays are stored as these can be easily converted into `diff` arrays but not the other way around.
        '''

        self.__get_variants()

        #Calculate differences
        self._nucleotides_full = self.__nucleotides()

        #These are only valid when a VCF has been applied to at least 1 of the genomes
        # self._indels_full = self.__indels()

        #Checking for the same genes, give a warning is the genes are different
        if self.genome1.genes.keys() != self.genome2.genes.keys():
            self.__raise_mutations_warning(self.genome1, self.genome2)

        self.update_view(self._view_method,'gene')

    def __snp_distance(self):
        '''Calculates the SNP distance between the two genomes

        Returns:
            int: The SNP distance between the two genomes
        '''
        # Ignores `z` and `x`
        return sum([1 for (base1, base2) in zip(self.genome1.nucleotide_sequence, self.genome2.nucleotide_sequence) if base1 != base2 and base1 in ['a','t','c','g'] and base2 in ['a','t','c','g']])

    def __get_variants(self):

        variants=[]
        indices=[]
        is_snp=[]
        is_het=[]
        is_null=[]
        is_indel=[]
        indel_length=[]
        indel_nucleotides=[]
        refs=[]
        alts=[]

        # first do the SNPs, HETs and NULLs
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence

        # for now we simply allow the ref and alt to be different e.g. can have a SNP on a NULL (x>t)
        for (r,idx,a) in zip(self.genome1.nucleotide_sequence[mask], self.genome1.nucleotide_index[mask], self.genome2.nucleotide_sequence[mask]):

            variants.append(str(idx)+r+'>'+a)
            refs.append(r)
            alts.append(a)
            indices.append(idx)
            is_indel.append(False)
            indel_length.append(0)
            indel_nucleotides.append(None)
            if a=='z':
                is_het.append(True)
                is_snp.append(False)
                is_null.append(False)
            elif a=='x':
                is_het.append(False)
                is_snp.append(False)
                is_null.append(True)
            elif a in ['a','t','c','g']:
                is_het.append(False)
                is_snp.append(True)
                is_null.append(False)

        # INDELs are trickier: we have to deal with the case where both genomes have an indel at the same position
        # if they are different we catch fire since that is too difficult to parse right now
        # if they are the same, then there is no difference
        # the other cases are the usual one of there being an indel on the RHS and then an indel on the LHS (but not the RHS)
        # for the latter we 'reverse' the indel e.g. 1000_ins_at - None becomes 1000_del_at
        # subtraction becomes "how do we go from the LHS to the RHS in a-b?"

        # for indels catch fire if both genomes have different indels at the same position
        assert numpy.sum(self.genome1.is_indel & self.genome2.is_indel & (self.genome1.indel_nucleotides!=self.genome2.indel_nucleotides))==0, 'both genomes have indels of different lengths at one or more of the same positions -- this cannot be easily resolved!'

        # the other case is where there is an identical indel at a position but this leads to a difference of zero!

        # if the indel is on the RHS, then it is unchanged
        mask = self.genome2.is_indel & (self.genome1.indel_nucleotides!=self.genome2.indel_nucleotides)

        for (idx,length,alt,r,a) in zip(self.genome1.nucleotide_index[mask],self.genome2.indel_length[mask], self.genome2.indel_nucleotides[mask],self.genome1.nucleotide_sequence, self.genome2.nucleotide_sequence):
            indices.append(idx)
            is_indel.append(True)
            is_het.append(False)
            is_snp.append(False)
            is_null.append(False)
            indel_length.append(length)
            indel_nucleotides.append(alt)
            if length>0:
                variants.append(str(idx)+'_ins_'+str(alt))
            else:
                variants.append(str(idx)+'_del_'+str(alt))

        # if the indel is on the LHS, then it is unchanged, then it needs 'reversing' since we are returning how to get to the RHS from the LHS hence we delete an insertion etc
        mask = self.genome1.is_indel & (self.genome1.indel_nucleotides!=self.genome2.indel_nucleotides)

        for (idx,length,alt,r,a) in zip(self.genome1.nucleotide_index[mask],self.genome1.indel_length[mask], self.genome1.indel_nucleotides[mask],self.genome1.nucleotide_sequence, self.genome2.nucleotide_sequence):
            indices.append(idx)
            is_indel.append(True)
            is_het.append(False)
            is_snp.append(False)
            is_null.append(False)
            length*=-1
            indel_length.append(length)
            indel_nucleotides.append(alt)
            if length>0:
                variants.append(str(idx)+'_ins_'+str(alt))
            else:
                variants.append(str(idx)+'_del_'+str(alt))

        self.variants=numpy.array(variants)
        self.nucleotide_index=numpy.array(indices)
        self.is_indel=numpy.array(is_indel)
        self.indel_length=numpy.array(indel_length)
        self.indel_nucleotides=numpy.array(indel_nucleotides)
        self.is_snp=numpy.array(is_snp)
        self.is_het=numpy.array(is_het)
        self.is_null=numpy.array(is_null)

    def __nucleotides(self):
        '''Calculate the difference in nucleotides
        Returns:
            numpy.array: Numpy array of tuples of (genome1_nucleotide, genome2_nucleotide)
        '''
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence
        return numpy.array(list(zip(self.genome1.nucleotide_sequence[mask], self.genome2.nucleotide_sequence[mask])))

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
            #This function is never called, but if in future it is, this line needs attention (self.indel_indicies doesn't exist)
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

class GeneDifference(Difference):
    '''Object to store the differences within genes. The view system is inherited from the Difference class.
    Set to `full` for arrays of tuple values where applicable.
    Set to `diff` for arrays of values from gene1 where the values vary. Default.
    This can be updated using the update_view function.

    Instance variables:
        gene1 (gumpy.Gene): Gene object 1
        gene2 (gumpy.Gene): Gene object 2
        nucleotides (numpy.array): Array of the nucleotides at which the genes differ. Format depends on the current view.
        mutations (numpy.array): Array of mutations in GARC between the two Gene objects.
        indel_indices (numpy.array): Array of nucleotide numbers where the indel lengths in the two genes differ.
        indels (numpy.array): Array of indel lengths where the indel lengths differ. Format depends on the current view.
        codons (numpy.array): Array of codons where the two Gene objects have different codons. Format depends on the current view.
        amino_acid_sequence (numpy.array): Array of amino acids where the two Gene objects have different amino acids. Format depends on the current view.
    Functions:
        amino_acid_variants(int) -> dict: Takes an amino acid index and returns a dictionary containing data from a vcf
                                            file (if applicable) for attributes such as calls, ref, and alt for all nucleotides
                                            within the codons for this amino acid index. If these genes do not code protien, returns {}
        nucleotide_variants(int) -> dict: Takes a nucleotide index and returns a dictionary containing data from a vcf
                                            file (if applicable) for attributes such as calls, ref, and alt at the given index
    Inherited functions:
        update_view(str) -> None: Used to change the viewing method for instance variables. Input values are either `diff` or `full`
    '''
    def __init__(self, gene1, gene2):
        '''
        Constructor. Takes in two gene objects and calculates the difference in a few areas such as nucleotides, codons, and amino acids.

        Args:
            gene1 (gumpy.Gene): Gene object 1
            gene2 (gumpy.Gene): Gene object 2
        '''

        assert isinstance(gene1, gumpy.Gene)
        assert isinstance(gene2, gumpy.Gene)

        if gene1.total_number_nucleotides != gene2.total_number_nucleotides:
            #The lengths of the genes are different so comparing them is meaningless
            warnings.warn("The two genes ("+gene1.name+" and "+gene2.name+") are different lengths, so comparision failed...", UserWarning)
            return None
        if gene1.name != gene2.name:
            warnings.warn("The two genes given have different names ("+gene1.name+", "+gene2.name+") but the same length, continuing...", UserWarning)
        if gene1.codes_protein != gene2.codes_protein:
            warnings.warn(f"The two genes given do not have the same protein coding for {gene1.name}: Gene1 = {gene1.codes_protein}, Gene2 = {gene2.codes_protein}, so comparison failed...", UserWarning)
            return None
        self.gene1 = gene1
        self.gene2 = gene2
        self.codes_protein=gene1.codes_protein
        self._view_method = "diff"

        self._nucleotides_full = self.__nucleotides()
        # self.mutations = self.__mutations()
        self.__get_mutations()

        self._codons_full = self.__codons()
        self._amino_acids_full = self.__amino_acids()

        self.update_view("diff", 'gene') #Use inherited method to set the view

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

    def __get_mutations(self):

        mutations=[]
        amino_acid_number=[]
        nucleotide_number=[]
        nucleotide_index=[]
        gene_position=[]
        is_cds=[]
        is_indel=[]
        is_promoter=[]
        indel_length=[]
        indel_nucleotides=[]
        ref_nucleotides=[]
        alt_nucleotides=[]
        is_snp=[]
        is_het=[]
        is_null=[]

        if self.codes_protein:

            mask=self.gene1.codons!=self.gene2.codons
            for (r,num,a,codon1,codon2) in zip( self.gene1.amino_acid_sequence[mask],\
                                                self.gene1.amino_acid_number[mask],\
                                                self.gene2.amino_acid_sequence[mask],\
                                                self.gene1.codons[mask],\
                                                self.gene2.codons[mask]):
                # ref.append(r)
                # alt.append(a)
                mutations.append(r+str(num)+a)
                if a == 'X':
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a=='Z':
                    is_null.append(False)
                    is_het.append(True)
                    is_snp.append(False)
                else:
                    is_null.append(False)
                    is_het.append(False)
                    is_snp.append(True)
                amino_acid_number.append(num)
                nucleotide_number.append(None)
                nucleotide_index.append(None)
                gene_position.append(num)
                is_cds.append(True)
                is_indel.append(False)
                is_promoter.append(False)
                indel_length.append(None)
                indel_nucleotides.append(None)
                ref_nucleotides.append(codon1)
                alt_nucleotides.append(codon2)

            mask=(self.gene1.nucleotide_sequence!=self.gene2.nucleotide_sequence) & (self.gene1.is_promoter)

            for (r,num,a,idx) in zip(   self.gene1.nucleotide_sequence[mask],\
                                    self.gene1.nucleotide_number[mask],\
                                    self.gene2.nucleotide_sequence[mask],\
                                    self.gene1.nucleotide_index[mask]):

                mutations.append(r+str(num)+a)
                # ref.append(r)
                # alt.append(a)
                amino_acid_number.append(None)
                nucleotide_number.append(num)
                nucleotide_index.append(idx)
                gene_position.append(num)
                is_cds.append(False)
                is_indel.append(False)
                is_promoter.append(True)
                indel_length.append(None)
                indel_nucleotides.append(None)
                ref_nucleotides.append(r)
                alt_nucleotides.append(a)
                if a == 'x':
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a=='z':
                    is_null.append(False)
                    is_het.append(True)
                    is_snp.append(False)
                else:
                    is_null.append(False)
                    is_het.append(False)
                    is_snp.append(True)

        else:

            mask=(self.gene1.nucleotide_sequence!=self.gene2.nucleotide_sequence)

            for (r,num,a,idx) in zip(   self.gene1.nucleotide_sequence[mask],\
                                    self.gene1.nucleotide_number[mask],\
                                    self.gene2.nucleotide_sequence[mask],\
                                    self.gene1.nucleotide_index[mask]):
                mutations.append(r+str(num)+a)
                # ref.append(r)
                # alt.append(a)
                amino_acid_number.append(None)
                nucleotide_number.append(num)
                nucleotide_index.append(idx)
                gene_position.append(num)
                if num>0:
                    is_cds.append(True)
                    is_promoter.append(False)
                else:
                    is_cds.append(False)
                    is_promoter.append(True)
                is_indel.append(False)
                indel_length.append(None)
                indel_nucleotides.append(None)
                ref_nucleotides.append(r)
                alt_nucleotides.append(a)
                if a == 'x':
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a=='z':
                    is_null.append(False)
                    is_het.append(True)
                    is_snp.append(False)
                else:
                    is_null.append(False)
                    is_het.append(False)
                    is_snp.append(True)

        # now let's do indels
        assert numpy.sum((self.gene1.is_indel & self.gene2.is_indel) & (self.gene1.indel_nucleotides!=self.gene2.indel_nucleotides))==0, 'both genes have different indels at one or more of the same positions -- this cannot be easily be resolved!'

        mask=self.gene2.is_indel & (self.gene1.indel_nucleotides!=self.gene2.indel_nucleotides)

        for (num,length,alt,idx) in zip(    self.gene1.nucleotide_number[mask],\
                                            self.gene2.indel_length[mask],\
                                            self.gene2.indel_nucleotides[mask],\
                                            self.gene1.nucleotide_index[mask]):
            # ref.append(None)
            # alt.append(None)
            amino_acid_number.append(None)
            nucleotide_number.append(num)
            nucleotide_index.append(idx)
            gene_position.append(num)
            if length>0:
                mutations.append(str(num)+'_ins_'+str(alt))
            else:
                mutations.append(str(num)+'_del_'+str(alt))
            if num>0:
                is_cds.append(True)
                is_promoter.append(False)
            else:
                is_cds.append(False)
                is_promoter.append(True)
            is_indel.append(True)
            indel_length.append(length)
            indel_nucleotides.append(alt)
            ref_nucleotides.append(None)
            alt_nucleotides.append(None)
            is_null.append(False)
            is_het.append(False)
            is_snp.append(False)

        mask=self.gene1.is_indel & (self.gene1.indel_nucleotides!=self.gene2.indel_nucleotides)

        for (num,length,alt,idx) in zip(    self.gene1.nucleotide_number[mask],\
                                        self.gene1.indel_length[mask],\
                                        self.gene1.indel_nucleotides[mask],\
                                        self.gene1.nucleotide_index[mask]):
            # ref.append(None)
            # alt.append(None)
            amino_acid_number.append(None)
            nucleotide_number.append(num)
            nucleotide_index.append(idx)
            gene_position.append(num)
            length*=-1
            if length>0:
                mutations.append(str(num)+'_ins_'+str(alt))
            else:
                mutations.append(str(num)+'_del_'+str(alt))
            if num>0:
                is_cds.append(True)
                is_promoter.append(False)
            else:
                is_cds.append(False)
                is_promoter.append(True)
            is_indel.append(True)
            indel_length.append(length)
            indel_nucleotides.append(alt)
            ref_nucleotides.append(None)
            alt_nucleotides.append(None)
            is_null.append(False)
            is_het.append(False)
            is_snp.append(False)

        self.mutations=numpy.array(mutations)
        self.amino_acid_number=numpy.array(amino_acid_number)
        self.nucleotide_number=numpy.array(nucleotide_number)
        self.nucleotide_index=numpy.array(nucleotide_index)
        self.gene_position=numpy.array(gene_position)
        self.is_cds=numpy.array(is_cds)
        self.is_promoter=numpy.array(is_promoter)
        self.is_indel=numpy.array(is_indel)
        self.indel_length=numpy.array(indel_length)
        self.indel_nucleotides=numpy.array(indel_nucleotides)
        self.ref_nucleotides=numpy.array(ref_nucleotides)
        self.alt_nucleotides=numpy.array(alt_nucleotides)
        self.is_snp=numpy.array(is_snp)
        self.is_het=numpy.array(is_het)
        self.is_null=numpy.array(is_null)

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
            if codon1 != codon2:
                aa_diff.append((aa1, aa2))
        #Mutations should be in order of aa->indel/promoter, so pad to match 
        fixed = aa_diff
        #Promoters
        to_add = sum(1 for p in numpy.logical_or(self.is_promoter, self.is_indel) if p)
        for _ in range(to_add):
            fixed.append((None, None))
        return numpy.array(fixed)



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

def collapse_inner_dict(dict_):
    '''Takes a dict which also contains 1 or more internal dictionaries
    Converts to a single 1D dictionary, ignoring key clashes

    Args:
        dict_ (dict): Dictionary to collapse
    Returns:
        dict: Collapsed dict
    '''
    fixed = {}
    for key in dict_.keys():
        if isinstance(dict_[key], dict):
            #Dict so unpack
            fixed = {**fixed, **dict_[key]}
        else:
            fixed[key] = dict_[key]
    return fixed
