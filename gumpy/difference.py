"""
Used to find differences between Genomes and Genes, as well as the impact a VCF has on 
a genome.

Abstract classes:
    * Difference - Abstract class providing the ability to change data views

Classes:
    * GenomeDifference
    * GeneDifference

Exceptions:
    * FailedComparison

Functions:
    convert_nucleotides_codons(numpy.array) -> numpy.array: Converts an array of 
        nucleotides to an array of codons.
    setup_codon_aa_dict() -> dict: Returns a dictionary mapping codon->amino_acid
"""

import warnings
from abc import ABC  # Python library for abstract classes
from typing import Dict, List, Tuple
import numpy
from tqdm import tqdm
import gumpy


class FailedComparison(Exception):
    """Exception to be raised in cases of a failed comparion"""

    def __init__(self, message: str):
        """Constructor

        Args:
            message (str): Error message
        """
        self.message = message
        super().__init__(self.message)


class Difference(ABC):
    """
    Abstract class used to provide the ability to switch between views.
    Inherited by GenomeDifference and GeneDifference.

    This should not be instantiated.
    """

    # Give some default values to appease the linter
    _nucleotides_full = numpy.array([])
    codes_protein = None
    _indels_full = numpy.array([])

    def update_view(self, method: str):
        """Update the viewing method. Can either be `diff` or `full`:

        * `diff`: Where applicable, variables return arrays of values from object1
                where they are not equal to object 2
        * `full`: Where applicable, variables return arrays of tuples showing
                (object1_val, object2_val) where values are not equal between objects.

        Args:
            method (str): Name of the viewing method. Must be within `['diff', 'full']`
        """
        assert isinstance(method, str), (
            "Invalid method " + str(method) + " of type " + str(type(method))
        )
        assert method in ["diff", "full"], "Invalid method: " + method

        if method == "full":
            self.nucleotides = self._nucleotides_full

        if method == "diff":
            # Convert the full arrays into diff arrays
            self.nucleotides = self.__full_to_diff(self._nucleotides_full)
        self._view_method = method

    def __full_to_diff(self, array: numpy.ndarray) -> numpy.ndarray:
        """Convert an array from a full view to a diff view

        Args:
            array (numpy.ndarray): Array of tuples of values
        Returns:
            numpy.ndarray: Array of values from object1
        """
        if len(array) > 0:
            return numpy.array([item2 for (item1, item2) in array], dtype=object)
        else:
            return numpy.array([])


class GenomeDifference(Difference):
    """
    GenomeDifference object captures the difference between two genomes. Other than
    `snp_distance`, all public instance variables should be arrays
    which align so `x[N] <--> y[N]`

    * The difference can be viewed in one of two ways:
        * `diff`: Arrays of the values from genome1 where there are different values
            in genome2, or values which exist in genome1 but not genome2 depending on
            which is more appropriate (closer to classical subtraction).
            This is the default view.
        * `full`: Arrays of tuples consisting of (genome1_val, genome2_val) -
            more information but less like classical subtraction
            and more difficult to wrangle meaningful information from.
        * This option can be set by using the update_view() method

    * Instance variables:
        * snp_distance (int): SNP distance between the two genomes
        * variants (numpy.array): Variants in GARC
        * nucleotide_index (numpy.array): Genome index (1 based) of this variant
        * indel_length (numpy.array): Length of the indel of this variant (None is
            not an indel)
        * indel_nucleotides (numpy.array): Nucleotides of this indel of this variant
            (None if not an indel)
        * vcf_evidence (numpy.array): VCF evidence to support this call
        * gene_name (numpy.array): Name of the gene this variant affects
            (None if intergene)
        * gene_position (numpy.array): Position of this variant within this gene.
            This refers to codon number if applicable (None if intergene)
        * codon_idx (numpy.array): 0-based index of the nucleotide within the codon
            this affects if applicable. None otherwise
    * Functions:
        * variants(int) -> dict: Takes a genome index and returns a dictionary mapping
            field->(genome1_val, genome2_val) for all fields of a vcf file
            (if applicable)
        * gene_differences() -> [GeneDifference]: Returns a list of GeneDifference
            objects
        * minor_populations() -> [str]: Returns a list of minor population mutations
            in GARC
    * Inherited functions:
        * update_view(str) -> None: Used to change the viewing method for instance
            variables. Input values are either `diff` or `full`
    """

    def __init__(self, genome1, genome2):
        """
        Constructor for the GenomeDifference object.

        Called implictly when `genome1.difference(genome2)` is performed.
        Args:
            genome1 (gumpy.Genome): A Genome object to compare against
            genome2 (gumpy.Genome): The other Genome object
        """

        # insist that both must be Genome objects
        assert isinstance(genome1, gumpy.Genome)
        assert isinstance(genome2, gumpy.Genome)

        # Checking for the same genes, give an error if the genes are different
        if genome1.genes.keys() != genome2.genes.keys():
            genes_in_1_not_2 = set(genome1.genes.keys()).difference(
                genome2.genes.keys()
            )
            genes_in_2_not_1 = set(genome2.genes.keys()).difference(
                genome1.genes.keys()
            )
            message = "The two genomes have different genes."
            if len(genes_in_1_not_2) > 0:
                message += (
                    f" Genome 1 has {len(genes_in_1_not_2)} more "
                    f"genes: {genes_in_1_not_2}."
                )
            if len(genes_in_2_not_1) > 0:
                message += (
                    f" Genome2 has {len(genes_in_2_not_1)} more "
                    f"genes: {genes_in_2_not_1}."
                )
            raise FailedComparison(message)

        self.genome1 = genome1
        self.genome2 = genome2
        self._view_method = "diff"

        # Calculate SNPs
        self.snp_distance = self.__snp_distance()
        """
        Where applicable, the `full` difference arrays are stored as these can be 
        easily converted into `diff` arrays but not the other way around.
        """

        self.__get_variants()

        # Calculate differences
        self._nucleotides_full = self.__nucleotides()

        self.__assign_vcf_evidence()

        self.update_view(self._view_method)

    def get_gene_pos(
        self, gene: str, idx: int, variant: str, start: int | None = None
    ) -> int | None:
        """Find the gene position of a given nucleotide index.
        This is considerably faster than building a whole stacked_gene_pos array
            (takes ~4.5mins for tb)

        Args:
            gene (str): Name of the gene to search
            idx (int): Nucleotide index we want the gene position of
            variant (str): Variant we're looking for in GARC
            start (int): Start position. Defaults to None

        Returns:
            int | None: Gene position of this nucleotide index
                        (or None if it lies outside of the gene)
        """
        stacked_gene_mask = self.genome1.stacked_gene_name == gene
        nc_idx = self.genome1.stacked_nucleotide_index[stacked_gene_mask]
        nc_nums = self.genome1.stacked_nucleotide_number[stacked_gene_mask]

        # Get the gene's nucleotide number
        nc_num = nc_nums[nc_idx == idx][0]

        # If coding SNP, pull out the codon number instead
        if (
            self.genome2.genes[gene]["codes_protein"]
            and nc_num > 0
            and "ins" not in variant
            and "del" not in variant
        ):
            # Use floor division (//)
            codon_number = (nc_num + 2) // 3
            return codon_number
        elif "ins" in variant:
            # Insertions need a little nudge in revcomp because of `ins at, del after`
            if self.genome2.genes[gene]["reverse_complement"]:
                nc_num -= 1
        elif "del" in variant:
            pos, t, bases = variant.split("_")
            if start is not None and idx != start:
                # This is not the start of the deletion, so truncate the bases as
                #   appropriate
                dels = len(bases) - (idx - start)
            else:
                dels = len(bases)
            # Deletions need even more nudging in revcomp because the entire deletion
            #   is reversed so starts at the end
            if self.genome2.genes[gene]["reverse_complement"]:
                nc_num = nc_num - dels + 1

                # Edge case of deletion starting in revcomp gene and extending
                # past gene start, so return None
                if (
                    self.genome2.stacked_nucleotide_number[
                        self.genome2.stacked_gene_name == gene
                    ][-1]
                    > nc_num
                ):
                    return None

        return nc_num

    def _get_vcf_idx(self, vcf_row: Dict) -> int | None:
        """Given a vcf row, figure out which alt it refers to. Should be available in
            the GT field

        Args:
            vcf_row (dict): Internal representation of the VCF row for this variant

        Returns:
            int | None: 0-based index of the alts this refers to
        """
        # Use `GT` field to figure out which one this call is
        # (This will be **much** messier for minor populations)
        gt = vcf_row["GT"]
        if gt[0] == gt[1]:
            # Actual call so return the 1-based idx
            return gt[0]
        else:
            # Het call so return None as it involves >1
            return None

    def __assign_vcf_evidence(self) -> None:
        """Once we have pulled out all of the variants, find the VCF evidence
        (if existing) for each.
        """
        evidences: List = []
        indices: List = []
        for i, idx in enumerate(self.nucleotide_index):
            evidence1 = self.genome1.vcf_evidence.get(idx)
            evidence2 = self.genome2.vcf_evidence.get(idx)

            # As we add the vcf evidence for every deleted base
            # This should only be reported for the actual start
            if (
                self.genome1.is_deleted[idx - 1]
                and self.genome1.indel_length[idx - 1] >= 0
            ):
                evidence1 = None
            if (
                self.genome2.is_deleted[idx - 1]
                and self.genome2.indel_length[idx - 1] >= 0
            ):
                evidence2 = None

            if evidence1 is not None and evidence2 is not None:
                # We have a collision. For now, just concat FIXME
                evidences.append([evidence1, evidence2])
                indices.append(
                    [self._get_vcf_idx(evidence1), self._get_vcf_idx(evidence2)]
                )
            elif evidence1 is not None:
                evidences.append(evidence1)
                indices.append(self._get_vcf_idx(evidence1))
            elif evidence2 is not None:
                evidences.append(evidence2)
                indices.append(self._get_vcf_idx(evidence2))
            else:
                evidences.append(None)
                indices.append(None)

        self.vcf_evidences = evidences
        self.vcf_idx = indices

    def minor_populations(self, interpretation: str = "reads") -> List[str]:
        """Get the minor population mutations in GARC

        Args:
            interpretation (str, optional): How to report minor population.
                'reads' reports number of reads. 'percentage' reports fractional
                read support. Defaults to 'reads'.

        Returns:
            List[str]: List of mutations in GARC
        """
        return self.genome2.minority_populations_GARC(
            interpretation=interpretation, reference=self.genome1
        )

    def __snp_distance(self) -> int:
        """Calculates the SNP distance between the two genomes

        Returns:
            int: The SNP distance between the two genomes
        """
        # Ignores `z` and `x`
        return sum(
            [
                1
                for (base1, base2) in zip(
                    self.genome1.nucleotide_sequence, self.genome2.nucleotide_sequence
                )
                if base1 != base2
                and base1 in ["a", "t", "c", "g"]
                and base2 in ["a", "t", "c", "g"]
            ]
        )

    def __get_variants(self):
        """Get the variants between the two genomes, populating the internal arrays"""

        variants = []
        indices = []
        is_snp = []
        is_het = []
        is_null = []
        is_indel = []
        indel_length = []
        indel_nucleotides = []
        refs = []
        alts = []
        gene_name = []
        gene_pos = []
        codon_idx = []

        # first do the SNPs, HETs and NULLs
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence

        # for now we simply allow the ref and alt to be different e.g. can have a
        #   SNP on a NULL (x>t)
        for r, idx, a in tqdm(
            zip(
                self.genome1.nucleotide_sequence[mask],
                self.genome1.nucleotide_index[mask],
                self.genome2.nucleotide_sequence[mask],
            ),
            total=len(self.genome1.nucleotide_sequence[mask]),
        ):
            variants.append(str(idx) + r + ">" + a)
            refs.append(r)
            alts.append(a)
            indices.append(idx)
            is_indel.append(False)
            indel_length.append(0)
            indel_nucleotides.append(None)
            if a == "z":
                is_het.append(True)
                is_snp.append(False)
                is_null.append(False)
            elif a == "x":
                is_het.append(False)
                is_snp.append(False)
                is_null.append(True)
            elif a in ["a", "t", "c", "g"]:
                is_het.append(False)
                is_snp.append(True)
                is_null.append(False)

            # Find the genes at this pos
            genes = sorted(
                list(
                    set(
                        self.genome1.stacked_gene_name[
                            self.genome1.stacked_nucleotide_index == idx
                        ]
                    )
                )
            )
            if len(genes) > 1:
                # If we have genes, we need to duplicate some info
                first = True
                for gene in genes:
                    if gene == "":
                        # If we have genes, we don't care about this one
                        continue
                    gene_name.append(gene)
                    gene_pos.append(self.get_gene_pos(gene, idx, variants[-1]))
                    if (
                        self.genome2.genes[gene]["codes_protein"]
                        and gene_pos[-1] is not None
                        and gene_pos[-1] > 0
                    ):
                        # Get codon idx
                        nc_idx = self.genome1.stacked_nucleotide_index[
                            self.genome1.stacked_gene_name == gene
                        ]
                        nc_num = self.genome1.stacked_nucleotide_number[
                            self.genome1.stacked_gene_name == gene
                        ]
                        codon_idx.append((nc_num[nc_idx == idx][0] - 1) % 3)
                    else:
                        codon_idx.append(None)

                    # If this isn't the first one, we need to duplicate the row
                    if first:
                        first = False
                    else:
                        variants.append(variants[-1])
                        refs.append(refs[-1])
                        alts.append(alts[-1])
                        indices.append(indices[-1])
                        is_indel.append(is_indel[-1])
                        indel_length.append(indel_length[-1])
                        indel_nucleotides.append(indel_nucleotides[-1])
                        is_het.append(is_het[-1])
                        is_snp.append(is_snp[-1])
                        is_null.append(is_null[-1])

            else:
                # We have 1 gene or none, so set to None if no gene is present
                gene = genes[0] if genes[0] != "" else None
                gene_name.append(gene)
                if gene is not None:
                    # Single gene, so pull out data
                    gene_pos.append(self.get_gene_pos(gene, idx, variants[-1]))

                    if (
                        self.genome2.genes[gene]["codes_protein"]
                        and gene_pos[-1] is not None
                        and gene_pos[-1] > 0
                    ):
                        # Get codon idx
                        nc_idx = self.genome1.stacked_nucleotide_index[
                            self.genome1.stacked_gene_name == gene
                        ]
                        nc_num = self.genome1.stacked_nucleotide_number[
                            self.genome1.stacked_gene_name == gene
                        ]
                        codon_idx.append((nc_num[nc_idx == idx][0] - 1) % 3)
                    else:
                        codon_idx.append(None)
                else:
                    gene_pos.append(None)
                    codon_idx.append(None)

        # INDELs are trickier: we have to deal with the case where both genomes have
        #   an indel at the same position
        # if they are different we catch fire since that is too difficult to parse
        #   right now
        # if they are the same, then there is no difference
        # the other cases are the usual one of there being an indel on the RHS and then
        #   an indel on the LHS (but not the RHS)
        # for the latter we 'reverse' the indel e.g. 1000_ins_at - None becomes 1000_del
        #   at
        # subtraction becomes "how do we go from the LHS to the RHS in a-b?"

        # for indels catch fire if both genomes have different indels at the same
        #   position
        assert (
            numpy.sum(
                self.genome1.is_indel
                & self.genome2.is_indel
                & (self.genome1.indel_nucleotides != self.genome2.indel_nucleotides)
            )
            == 0
        ), (
            "both genomes have indels of different lengths at one or more of the "
            "same positions -- this cannot be easily resolved!"
        )

        # the other case is where there is an identical indel at a position but
        #   this leads to a difference of zero!

        # if the indel is on the RHS, then it is unchanged
        mask = self.genome2.is_indel & (
            self.genome1.indel_nucleotides != self.genome2.indel_nucleotides
        )

        for idx, length, alt, r, a in zip(
            self.genome1.nucleotide_index[mask],
            self.genome2.indel_length[mask],
            self.genome2.indel_nucleotides[mask],
            self.genome1.nucleotide_sequence,
            self.genome2.nucleotide_sequence,
        ):
            indices.append(idx)
            is_indel.append(True)
            is_het.append(False)
            is_snp.append(False)
            is_null.append(False)
            refs.append(None)
            alts.append(None)
            indel_length.append(length)
            indel_nucleotides.append(alt)
            if length > 0:
                variants.append(str(idx) + "_ins_" + str(alt))
            else:
                variants.append(str(idx) + "_del_" + str(alt))

            # Find the genes at this pos
            positions = []
            seen_genes = set()
            if "ins" in variants[-1]:
                # Insertion, so only need to check a single index
                positions.append(idx)
                start = None
            else:
                # Deletion, so we need to check every deleted position for genes
                pos, _, bases = variants[-1].split("_")
                start = idx
                for i in range(len(bases)):
                    positions.append(idx + i)
            first = True
            for pos in positions:
                genes = sorted(
                    list(
                        set(
                            self.genome1.stacked_gene_name[
                                self.genome1.stacked_nucleotide_index == pos
                            ]
                        )
                    )
                )
                adding = False
                if len(genes) > 1:
                    multi = True
                    # If we have genes, we need to duplicate some info
                    for gene in genes:
                        if gene == "":
                            # If we have genes, we don't care about this one
                            continue

                        if gene in seen_genes:
                            # If we've already seen it, skip
                            continue
                        else:
                            # Mark as seen now we have seen it
                            seen_genes.add(gene)
                            adding = True
                        gene_name.append(gene)
                        gene_pos.append(
                            self.get_gene_pos(gene, pos, variants[-1], start=start)
                        )
                        if (
                            self.genome2.genes[gene]["codes_protein"]
                            and gene_pos[-1] is not None
                            and gene_pos[-1] > 0
                        ):
                            # Get codon pos
                            nc_idx = self.genome1.stacked_nucleotide_index[
                                self.genome1.stacked_gene_name == gene
                            ]
                            nc_num = self.genome1.stacked_nucleotide_number[
                                self.genome1.stacked_gene_name == gene
                            ]
                            codon_idx.append((nc_num[nc_idx == pos][0] - 1) % 3)
                        else:
                            codon_idx.append(None)

                        # If this isn't the first one, we need to duplicate the row
                        if first:
                            first = False
                        else:
                            variants.append(variants[-1])
                            refs.append(refs[-1])
                            alts.append(alt[-1])
                            indices.append(indices[-1])
                            is_indel.append(is_indel[-1])
                            indel_length.append(indel_length[-1])
                            indel_nucleotides.append(indel_nucleotides[-1])
                            is_het.append(is_het[-1])
                            is_snp.append(is_snp[-1])
                            is_null.append(is_null[-1])
                else:
                    # We have 1 gene or none, so set to None if no gene is present
                    gene = genes[0] if genes[0] != "" else None
                    if gene in seen_genes:
                        # If we've already seen it, or it's intragene so skip
                        continue
                    else:
                        # Mark as seen now we have seen it
                        seen_genes.add(gene)
                        adding = True
                    if gene is not None:
                        # Single gene, so pull out data
                        gene_name.append(gene)
                        gene_pos.append(
                            self.get_gene_pos(gene, pos, variants[-1], start=start)
                        )

                        if (
                            self.genome2.genes[gene]["codes_protein"]
                            and gene_pos[-1] is not None
                            and gene_pos[-1] > 0
                        ):
                            # Get codon pos
                            nc_idx = self.genome1.stacked_nucleotide_index[
                                self.genome1.stacked_gene_name == gene
                            ]
                            nc_num = self.genome1.stacked_nucleotide_number[
                                self.genome1.stacked_gene_name == gene
                            ]
                            codon_idx.append((nc_num[nc_idx == pos][0] - 1) % 3)
                        else:
                            codon_idx.append(None)
                    else:
                        # We don't care about intragene values if its not the first item
                        if pos != idx:
                            continue
                        gene_name.append(None)
                        gene_pos.append(None)
                        codon_idx.append(None)
                # If this isn't the first one, we need to duplicate the row
                if first:
                    first = False
                elif adding and not multi:
                    variants.append(variants[-1])
                    refs.append(refs[-1])
                    alts.append(alt[-1])
                    indices.append(indices[-1])
                    is_indel.append(is_indel[-1])
                    indel_length.append(indel_length[-1])
                    indel_nucleotides.append(indel_nucleotides[-1])
                    is_het.append(is_het[-1])
                    is_snp.append(is_snp[-1])
                    is_null.append(is_null[-1])

        # if the indel is on the LHS, then it is unchanged, then it needs 'reversing'
        #   since we are returning how to get to the RHS from the LHS hence we
        #   delete an insertion etc
        mask = self.genome1.is_indel & (
            self.genome1.indel_nucleotides != self.genome2.indel_nucleotides
        )

        for idx, length, alt, r, a in zip(
            self.genome1.nucleotide_index[mask],
            self.genome1.indel_length[mask],
            self.genome1.indel_nucleotides[mask],
            self.genome1.nucleotide_sequence,
            self.genome2.nucleotide_sequence,
        ):
            indices.append(idx)
            is_indel.append(True)
            is_het.append(False)
            is_snp.append(False)
            is_null.append(False)
            refs.append(None)
            alts.append(None)
            length *= -1
            indel_length.append(length)
            indel_nucleotides.append(alt)
            if length > 0:
                variants.append(str(idx) + "_ins_" + str(alt))
            else:
                variants.append(str(idx) + "_del_" + str(alt))

            # Find the genes at this pos
            positions = []
            seen_genes = set()
            if "ins" in variants[-1]:
                # Insertion, so only need to check a single index
                positions.append(idx)
                start = None
            else:
                # Deletion, so we need to check every deleted position for genes
                pos, _, bases = variants[-1].split("_")
                start = idx
                for i in range(len(bases)):
                    positions.append(idx + i)
            first = True
            for pos in positions:
                genes = sorted(
                    list(
                        set(
                            self.genome1.stacked_gene_name[
                                self.genome1.stacked_nucleotide_index == pos
                            ]
                        )
                    )
                )
                adding = False
                multi = False
                if len(genes) > 1:
                    multi = True
                    # If we have genes, we need to duplicate some info
                    for gene in genes:
                        if gene == "":
                            # If we have genes, we don't care about this one
                            continue

                        if gene in seen_genes:
                            # If we've already seen it, skip
                            continue
                        else:
                            # Mark as seen now we have seen it
                            seen_genes.add(gene)
                            adding = True

                        gene_name.append(gene)
                        gene_pos.append(
                            self.get_gene_pos(gene, pos, variants[-1], start=start)
                        )
                        if (
                            self.genome2.genes[gene]["codes_protein"]
                            and gene_pos[-1] is not None
                            and gene_pos[-1] > 0
                        ):
                            # Get codon pos
                            nc_idx = self.genome1.stacked_nucleotide_index[
                                self.genome1.stacked_gene_name == gene
                            ]
                            nc_num = self.genome1.stacked_nucleotide_number[
                                self.genome1.stacked_gene_name == gene
                            ]
                            codon_idx.append((nc_num[nc_idx == pos][0] - 1) % 3)
                        else:
                            codon_idx.append(None)
                        # If this isn't the first one, we need to duplicate the row
                        if first:
                            first = False
                        else:
                            variants.append(variants[-1])
                            refs.append(refs[-1])
                            alts.append(alt[-1])
                            indices.append(indices[-1])
                            is_indel.append(is_indel[-1])
                            indel_length.append(indel_length[-1])
                            indel_nucleotides.append(indel_nucleotides[-1])
                            is_het.append(is_het[-1])
                            is_snp.append(is_snp[-1])
                            is_null.append(is_null[-1])

                else:
                    # We have 1 gene or none, so set to None if no gene is present
                    gene = genes[0] if genes[0] != "" else None
                    if gene in seen_genes:
                        # If we've already seen it, or it's intragene so skip
                        continue
                    else:
                        # Mark as seen now we have seen it
                        seen_genes.add(gene)
                        adding = True
                    if gene is not None:
                        # Single gene, so pull out data
                        gene_name.append(gene)
                        gene_pos.append(
                            self.get_gene_pos(gene, pos, variants[-1], start=start)
                        )

                        if (
                            self.genome2.genes[gene]["codes_protein"]
                            and gene_pos[-1] is not None
                            and gene_pos[-1] > 0
                        ):
                            # Get codon pos
                            nc_idx = self.genome1.stacked_nucleotide_index[
                                self.genome1.stacked_gene_name == gene
                            ]
                            nc_num = self.genome1.stacked_nucleotide_number[
                                self.genome1.stacked_gene_name == gene
                            ]
                            codon_idx.append((nc_num[nc_idx == pos][0] - 1) % 3)
                        else:
                            codon_idx.append(None)
                    else:
                        # We don't care about intragene values if its not the first item
                        if pos != idx:
                            continue
                        gene_name.append(None)
                        gene_pos.append(None)
                        codon_idx.append(None)

                # If this isn't the first one, we need to duplicate the row
                if first:
                    first = False
                elif adding and not multi:
                    variants.append(variants[-1])
                    refs.append(refs[-1])
                    alts.append(alt[-1])
                    indices.append(indices[-1])
                    is_indel.append(is_indel[-1])
                    indel_length.append(indel_length[-1])
                    indel_nucleotides.append(indel_nucleotides[-1])
                    is_het.append(is_het[-1])
                    is_snp.append(is_snp[-1])
                    is_null.append(is_null[-1])

        self.variants = numpy.array(variants)
        self.nucleotide_index = numpy.array(indices)
        self.is_indel = numpy.array(is_indel)
        self.indel_length = numpy.array(indel_length)
        self.indel_nucleotides = numpy.array(indel_nucleotides)
        self.is_snp = numpy.array(is_snp)
        self.is_het = numpy.array(is_het)
        self.is_null = numpy.array(is_null)
        self.gene_name = gene_name
        self.gene_pos = gene_pos
        self.codon_idx = codon_idx

    def __nucleotides(self) -> numpy.ndarray:
        """Calculate the difference in nucleotides
        Returns:
            numpy.ndarray: Numpy array of tuples of (genome1_nucleotide,
                genome2_nucleotide)
        """
        mask = self.genome1.nucleotide_sequence != self.genome2.nucleotide_sequence
        return numpy.array(
            list(
                zip(
                    self.genome1.nucleotide_sequence[mask],
                    self.genome2.nucleotide_sequence[mask],
                )
            )
        )


class GeneDifference(Difference):
    """Object to store the differences within genes. The view system is inherited from
        the Difference class.

    * Set to `full` for arrays of tuple values where applicable.
    * Set to `diff` for arrays of values from gene1 where the values vary. Default.
    * This can be updated using the update_view function.

    * Instance variables:
        * gene1 (gumpy.Gene): Gene object 1
        * gene2 (gumpy.Gene): Gene object 2
        * nucleotides (numpy.array): Array of the nucleotides at which the genes differ.
            Format depends on the current view.
        * mutations (numpy.array): Array of mutations in GARC between the two Gene
            objects.
        * indels (numpy.array): Array of indel lengths where the indel lengths differ.
            Format depends on the current view.
        * codons (numpy.array): Array of codons where the two Gene objects have
            different codons. Format depends on the current view.
        * amino_acid_sequence (numpy.array): Array of amino acids where the two Gene
            objects have different amino acids. Format depends on the current view.
    * Functions:
        * amino_acid_variants(int) -> dict: Takes an amino acid index and returns a
            dictionary containing data from a vcf file (if applicable) for attributes
            such as calls, ref, and alt for all nucleotides within the codons for this
            amino acid index. If these genes do not code protien, returns {}
        * nucleotide_variants(int) -> dict: Takes a nucleotide index and returns a
            dictionary containing data from a vcf file (if applicable) for attributes
            such as calls, ref, and alt at the given index
        * minor_populations() -> [str]: Returns a list of minor population mutations
            in GARC
    * Inherited functions:
        * update_view(str) -> None: Used to change the viewing method for instance
            variables. Input values are either `diff` or `full`
    """

    def __init__(self, gene1, gene2):
        """
        Constructor. Takes in two gene objects and calculates the difference in a few
            areas such as nucleotides, codons, and amino acids.

        Args:
            gene1 (gumpy.Gene): Gene object 1
            gene2 (gumpy.Gene): Gene object 2
        """

        assert isinstance(gene1, gumpy.Gene)
        assert isinstance(gene2, gumpy.Gene)

        if gene1.total_number_nucleotides != gene2.total_number_nucleotides:
            # The lengths of the genes are different so comparing them is meaningless
            raise FailedComparison(
                "The two genes ("
                + gene1.name
                + " and "
                + gene2.name
                + ") are different lengths, so comparision failed..."
            )
        if gene1.name != gene2.name:
            warnings.warn(
                "The two genes given have different names ("
                + gene1.name
                + ", "
                + gene2.name
                + ") but the same length, continuing...",
                UserWarning,
            )
        if gene1.codes_protein != gene2.codes_protein:
            raise FailedComparison(
                f"The two genes given do not have the same protein coding "
                f"for {gene1.name}: Gene1 = {gene1.codes_protein}, "
                f"Gene2 = {gene2.codes_protein}, so comparison failed..."
            )
        self.gene1 = gene1
        self.gene2 = gene2
        self.codes_protein = gene1.codes_protein
        self._view_method = "diff"

        self._nucleotides_full = self.__nucleotides()
        # self.mutations = self.__mutations()
        self.__get_mutations()

        self._codons_full = self.__codons()
        self._amino_acids_full = self.__amino_acids()

        self.__large_deletions()
        self.__assign_vcf_evidence()

        self.update_view("diff")  # Use inherited method to set the view

    def __assign_vcf_evidence(self) -> None:
        """Once we have pulled out all of the mutations, find the VCF evidence
        (if existing) for each
        """
        evidences: List = []
        for idx in self.nucleotide_index:
            evidence1 = self.gene1.vcf_evidence.get(idx)
            evidence2 = self.gene2.vcf_evidence.get(idx)

            if evidence1 is not None and evidence2 is not None:
                # We have a collision. For now, just concat FIXME
                evidences.append([evidence1, evidence2])
            elif evidence1 is not None:
                evidences.append(evidence1)
            elif evidence2 is not None:
                evidences.append(evidence2)
            else:
                evidences.append(None)
        self.vcf_evidences = evidences

    def __large_deletions(self):
        """Check to see what proportion of the gene has been deleted. Report separately
            than other mutations if >=50% deleted

        Updates internal arrays to denote a new mutation if deletion above threshold.
        """
        # Find everywhere that we have different deletion within coding regions
        mask = self.gene1.is_deleted != self.gene2.is_deleted

        # Find out how much of the gene we have deleted
        total = sum(self.gene1.is_deleted[mask]) + sum(self.gene2.is_deleted[mask])

        if total > 0:
            # We have some deletions
            percentage = total / len(self.gene1.nucleotide_sequence)
            if percentage >= 0.5:
                # More than 50% deleted, so give a percentage
                self.mutations = numpy.append(
                    self.mutations, [f"del_{round(percentage, 2)}"]
                )
                # Pull out the start of the deletion for vcf evidence
                self.nucleotide_index = numpy.append(
                    self.nucleotide_index, [self.gene2.nucleotide_index[mask][0]]
                )

                self.amino_acid_number = numpy.append(self.amino_acid_number, [None])

                # Get the nucleotide number of the first base of the gene
                # Not predictable as promoters are variable length
                self.nucleotide_number = numpy.append(
                    self.nucleotide_number,
                    [
                        self.gene2.nucleotide_number[
                            self.gene1.is_deleted | self.gene2.is_deleted
                        ][0]
                    ],
                )

                self.gene_position = numpy.append(self.gene_position, [None])
                self.is_cds = numpy.append(self.is_cds, [True])
                self.is_promoter = numpy.append(self.is_promoter, [False])
                self.is_indel = numpy.append(self.is_indel, [True])
                self.indel_length = numpy.append(self.indel_length, [-1 * total])
                self.indel_nucleotides = numpy.append(
                    self.indel_nucleotides,
                    ["".join(self.gene1.nucleotide_sequence[mask])],
                )
                self.ref_nucleotides = numpy.append(self.ref_nucleotides, [None])
                self.alt_nucleotides = numpy.append(self.alt_nucleotides, [None])
                self.is_snp = numpy.append(self.is_snp, [False])
                self.is_het = numpy.append(self.is_het, [False])
                self.is_null = numpy.append(self.is_null, [False])
            else:
                # Check for a deletion at the start of the gene which is the result of
                #   an upstream del
                if self.gene1.reverse_complement:
                    # This is difficult as revcomp shifts dels upstream by the length
                    #   of the del
                    # So we have to find a del in the RHS of the gene which would be
                    #   the last item if not shifted
                    pos = -1
                    while pos * -1 < len(self.gene1.nucleotide_sequence) / 2:
                        if (
                            self.gene1.indel_length[pos] < 0
                            and not self.gene1.is_indel[pos]
                        ):
                            # Del here in gene1 so check for validity
                            if (
                                self.gene1.nucleotide_index[pos]
                                + self.gene1.indel_length[pos]
                                + 1
                                == self.gene1.nucleotide_index[-1]
                            ):
                                # Re-adusting for deletion, this matches
                                break
                        if (
                            self.gene2.indel_length[pos] < 0
                            and not self.gene2.is_indel[pos]
                        ):
                            # Del here in gene2 so check for validity
                            if (
                                self.gene2.nucleotide_index[pos]
                                + self.gene2.indel_length[pos]
                                + 1
                                == self.gene2.nucleotide_index[-1]
                            ):
                                # Re-adusting for deletion, this matches
                                break
                        pos -= 1
                else:
                    # If not revcomp, the del is the first item
                    pos = 0

                if (
                    self.gene1.indel_length[pos] < 0 and not self.gene1.is_indel[pos]
                ) and (
                    self.gene2.indel_length[pos] < 0 and not self.gene2.is_indel[pos]
                ):
                    # Both have a del, so no difference
                    return

                if self.gene1.indel_length[pos] < 0 and not self.gene1.is_indel[pos]:
                    deleted = self.gene1.indel_nucleotides[pos]
                elif self.gene2.indel_length[pos] < 0 and not self.gene2.is_indel[pos]:
                    deleted = self.gene2.indel_nucleotides[pos]
                else:
                    # Neither, so no difference either
                    return

                # Adjust idx to point to the  right VCF evidence
                if self.gene1.reverse_complement:
                    idx = self.gene1.nucleotide_index[pos + len(deleted) - 1]
                else:
                    idx = self.gene1.nucleotide_index[pos]

                pos = self.gene1.nucleotide_number[pos]

                self.mutations = numpy.append(self.mutations, [f"{pos}_del_{deleted}"])
                # Pull out the start of the deletion for vcf evidence
                self.nucleotide_index = numpy.append(self.nucleotide_index, [idx])
                self.amino_acid_number = numpy.append(self.amino_acid_number, [None])

                # Get the nucleotide number of the first base of the gene
                # Not predictable as promoters are variable length
                self.nucleotide_number = numpy.append(
                    self.nucleotide_number,
                    [
                        self.gene2.nucleotide_number[
                            self.gene1.is_deleted | self.gene2.is_deleted
                        ][0]
                    ],
                )

                self.gene_position = numpy.append(self.gene_position, [pos])
                self.is_cds = numpy.append(self.is_cds, [True])
                self.is_promoter = numpy.append(self.is_promoter, [False])
                self.is_indel = numpy.append(self.is_indel, [True])
                self.indel_length = numpy.append(self.indel_length, [len(deleted)])
                self.indel_nucleotides = numpy.append(self.indel_nucleotides, [deleted])
                self.ref_nucleotides = numpy.append(self.ref_nucleotides, [None])
                self.alt_nucleotides = numpy.append(self.alt_nucleotides, [None])
                self.is_snp = numpy.append(self.is_snp, [False])
                self.is_het = numpy.append(self.is_het, [False])
                self.is_null = numpy.append(self.is_null, [False])

    def minor_populations(self, interpretation: str = "reads") -> List[str]:
        """Get the minor population mutations in GARC

        Args:
            interpretation (str, optional): How to report minor population. 'reads'
                reports number of reads. 'percentage' reports fractional read support.
                Defaults to 'reads'.

        Returns:
            [str]: List of mutations in GARC
        """
        return self.gene2.minority_populations_GARC(
            interpretation=interpretation, reference=self.gene1
        )

    def __nucleotides(self) -> numpy.ndarray:
        """Find the differences in nucleotides

        Returns:
            numpy.ndarray: Array of tuples (gene1_nucleotide, gene2_nucleotide)
        """
        return numpy.array(
            [
                (n1, n2)
                for (n1, n2) in zip(
                    self.gene1.nucleotide_sequence, self.gene2.nucleotide_sequence
                )
                if n1 != n2
            ]
        )

    def __get_mutations(self):
        """Get the mutations between the two genes, populating internal arrays"""

        mutations = []
        amino_acid_number = []
        nucleotide_number = []
        nucleotide_index = []
        gene_position = []
        is_cds = []
        is_indel = []
        is_promoter = []
        indel_length = []
        indel_nucleotides = []
        ref_nucleotides = []
        alt_nucleotides = []
        is_snp = []
        is_het = []
        is_null = []
        variants = []

        if self.codes_protein:
            mask = self.gene1.codons != self.gene2.codons
            for r, num, a, codon1, codon2 in zip(
                self.gene1.amino_acid_sequence[mask],
                self.gene1.amino_acid_number[mask],
                self.gene2.amino_acid_sequence[mask],
                self.gene1.codons[mask],
                self.gene2.codons[mask],
            ):
                # ref.append(r)
                # alt.append(a)
                mutations.append(r + str(num) + a)
                if a == "X":
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a == "Z":
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

                # Reconstruct the variant which caused this mutation to allow joining
                #   mutation and variant tables
                codon_nucleotide_indices = self.gene1.nucleotide_index[
                    self.gene1.is_cds
                ][self.gene1.codon_number == num]
                v = ""
                for i, (r_, a_) in enumerate(zip(codon1, codon2)):
                    if r_ != a_:
                        v += str(codon_nucleotide_indices[i]) + r_ + ">" + a_ + "&"
                if v[-1] == "&":
                    v = v[:-1]
                variants.append(v)

                # If synonymous mutation, pull out nucleotide variants too
                # This lets us determine effects of mutations such as fabG1@L203L more
                #   precisely
                if r == a:
                    for i, (rn, an) in enumerate(zip(codon1, codon2)):
                        if rn != an:
                            mutations.append(rn + str((num - 1) * 3 + i + 1) + an)
                            is_null.append(an == "x")
                            is_het.append(an == "z")
                            is_snp.append(an not in ["x", "z"])
                            amino_acid_number.append(None)
                            nucleotide_number.append((num - 1) * 3 + i + 1)
                            nucleotide_index.append(
                                self.gene1.nucleotide_index[
                                    self.gene1.nucleotide_number
                                    == (num - 1) * 3 + i + 1
                                ][0]
                            )
                            gene_position.append((num - 1) * 3 + i + 1)
                            is_cds.append(True)
                            is_indel.append(False)
                            is_promoter.append(False)
                            indel_length.append(None)
                            indel_nucleotides.append(None)
                            ref_nucleotides.append(rn)
                            alt_nucleotides.append(an)

                            variants.append(
                                str(codon_nucleotide_indices[i]) + rn + ">" + an
                            )

            mask = (
                self.gene1.nucleotide_sequence != self.gene2.nucleotide_sequence
            ) & (self.gene1.is_promoter)

            for r, num, a, idx in zip(
                self.gene1.nucleotide_sequence[mask],
                self.gene1.nucleotide_number[mask],
                self.gene2.nucleotide_sequence[mask],
                self.gene1.nucleotide_index[mask],
            ):
                mutations.append(r + str(num) + a)
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
                if a == "x":
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a == "z":
                    is_null.append(False)
                    is_het.append(True)
                    is_snp.append(False)
                else:
                    is_null.append(False)
                    is_het.append(False)
                    is_snp.append(True)

                # Reconstruct the variant which caused this mutation to allow joining
                #   mutation and variant tables
                variants.append(str(idx) + r + ">" + a)

        else:
            mask = self.gene1.nucleotide_sequence != self.gene2.nucleotide_sequence

            for r, num, a, idx in zip(
                self.gene1.nucleotide_sequence[mask],
                self.gene1.nucleotide_number[mask],
                self.gene2.nucleotide_sequence[mask],
                self.gene1.nucleotide_index[mask],
            ):
                mutations.append(r + str(num) + a)
                # ref.append(r)
                # alt.append(a)
                amino_acid_number.append(None)
                nucleotide_number.append(num)
                nucleotide_index.append(idx)
                gene_position.append(num)
                if num > 0:
                    is_cds.append(False)
                    is_promoter.append(False)
                else:
                    is_cds.append(False)
                    is_promoter.append(True)
                is_indel.append(False)
                indel_length.append(None)
                indel_nucleotides.append(None)
                ref_nucleotides.append(r)
                alt_nucleotides.append(a)
                if a == "x":
                    is_null.append(True)
                    is_het.append(False)
                    is_snp.append(False)
                elif a == "z":
                    is_null.append(False)
                    is_het.append(True)
                    is_snp.append(False)
                else:
                    is_null.append(False)
                    is_het.append(False)
                    is_snp.append(True)

                # Reconstruct the variant which caused this mutation to allow joining
                #   mutation and variant tables
                variants.append(str(idx) + r + ">" + a)

        # now let's do indels
        assert (
            numpy.sum(
                (self.gene1.is_indel & self.gene2.is_indel)
                & (self.gene1.indel_nucleotides != self.gene2.indel_nucleotides)
            )
            == 0
        ), (
            "both genes have different indels at one or more of the same positions -- "
            "this cannot be easily be resolved!"
        )

        mask = self.gene2.is_indel & (
            self.gene1.indel_nucleotides != self.gene2.indel_nucleotides
        )

        for num, length, alt, idx in zip(
            self.gene1.nucleotide_number[mask],
            self.gene2.indel_length[mask],
            self.gene2.indel_nucleotides[mask],
            self.gene2.indel_index[mask],
        ):
            # ref.append(None)
            # alt.append(None)
            amino_acid_number.append(None)
            nucleotide_number.append(num)
            nucleotide_index.append(idx)
            gene_position.append(num)
            if length > 0:
                mutations.append(str(num) + "_ins_" + str(alt))
            else:
                mutations.append(str(num) + "_del_" + str(alt))
            if num > 0:
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

            # Reconstruct the variant which caused this mutation to allow joining
            #   mutation and variant tables
            m = mutations[-1].split("_")[1:]
            if self.gene2.reverse_complement:
                # These are slightly different as revcomp need changes
                variants.append(
                    "_".join([self.gene2.revcomp_indel_nc_index[str(num)]] + m)
                )
            else:
                variants.append(
                    "_".join(
                        [
                            str(i)
                            for i in self.gene2.nucleotide_index[
                                self.gene2.nucleotide_number == num
                            ]
                        ]
                        + m
                    )
                )

        mask = self.gene1.is_indel & (
            self.gene1.indel_nucleotides != self.gene2.indel_nucleotides
        )

        for num, length, alt, idx in zip(
            self.gene1.nucleotide_number[mask],
            self.gene1.indel_length[mask],
            self.gene1.indel_nucleotides[mask],
            self.gene1.indel_index[mask],
        ):
            # ref.append(None)
            # alt.append(None)
            amino_acid_number.append(None)
            nucleotide_number.append(num)
            nucleotide_index.append(idx)
            gene_position.append(num)
            length *= -1
            if length > 0:
                mutations.append(str(num) + "_ins_" + str(alt))
            else:
                mutations.append(str(num) + "_del_" + str(alt))
            if num > 0:
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

            # Reconstruct the variant which caused this mutation to allow joining
            #   mutation and variant tables
            m = mutations[-1].split("_")[1:]
            if self.gene1.reverse_complement:
                # These are slightly different as revcomp need changes
                variants.append(
                    "_".join([self.gene1.revcomp_indel_nc_index[str(num)]] + m)
                )
            else:
                variants.append(
                    "_".join(
                        [
                            str(i)
                            for i in self.gene1.nucleotide_index[
                                self.gene1.nucleotide_number == num
                            ]
                        ]
                        + m
                    )
                )
        self.mutations = numpy.array(mutations)
        self.amino_acid_number = numpy.array(amino_acid_number)
        self.nucleotide_number = numpy.array(nucleotide_number)
        self.nucleotide_index = numpy.array(nucleotide_index)
        self.gene_position = numpy.array(gene_position)
        self.is_cds = numpy.array(is_cds)
        self.is_promoter = numpy.array(is_promoter)
        self.is_indel = numpy.array(is_indel)
        self.indel_length = numpy.array(indel_length)
        self.indel_nucleotides = numpy.array(indel_nucleotides)
        self.ref_nucleotides = numpy.array(ref_nucleotides)
        self.alt_nucleotides = numpy.array(alt_nucleotides)
        self.is_snp = numpy.array(is_snp)
        self.is_het = numpy.array(is_het)
        self.is_null = numpy.array(is_null)
        self.variants = variants

    def __codons(self) -> numpy.ndarray:
        """Find the codon positions which are different within the genes
            (within codon regions)

        Returns:
            numpy.ndarray: Array of codons which differ of the form
                [(gene1_codon, gene2_codon)]
        """
        codons1 = convert_nucleotides_codons(
            self.gene1.nucleotide_sequence[self.gene1.is_cds]
        )
        codons2 = convert_nucleotides_codons(
            self.gene2.nucleotide_sequence[self.gene2.is_cds]
        )
        mask = codons1 != codons2
        codons1 = codons1[mask]
        codons2 = codons2[mask]
        return numpy.array(list(zip(codons1, codons2)))

    def __amino_acids(self) -> numpy.ndarray:
        """Calculate the difference in amino acid sequences (within codon regions)

        Returns:
            numpy.ndarray: Array of tuples showing (amino_acid_1, amino_acid_2)
        """
        codon_to_amino_acid = setup_codon_aa_dict()

        aa_diff: List[Tuple[str | None, str | None]] = []
        for codon1, codon2 in self._codons_full:
            aa1 = codon_to_amino_acid[codon1]
            aa2 = codon_to_amino_acid[codon2]
            if codon1 != codon2:
                aa_diff.append((aa1, aa2))
                # Adding nucleotide variants for synon
                if aa1 == aa2:
                    for c1, c2 in zip(codon1, codon2):
                        if c1 != c2:
                            aa_diff.append((None, None))
        # Mutations should be in order of aa->indel/promoter, so pad to match
        fixed = aa_diff
        # Promoters
        to_add = sum(1 for p in numpy.logical_or(self.is_promoter, self.is_indel) if p)
        for _ in range(to_add):
            fixed.append((None, None))
        return numpy.array(fixed)


"""
Helper functions not specific to a class
"""


def convert_nucleotides_codons(nucleotides: numpy.ndarray) -> numpy.ndarray:
    """Helper function to convert an array of nucleotides into an array of codons

    Args:
        nucleotides (numpy.ndarray): Array of nucleotides

    Returns:
        numpy.ndarray: Array of codons
    """
    codons = []
    c = ""
    for index in range(1, len(nucleotides) + 1):
        c += nucleotides[index - 1]
        if index % 3 == 0:
            # There have been 3 bases seen so add the codon
            codons.append(c)
            c = ""
    return numpy.array(codons)


def setup_codon_aa_dict() -> Dict[str, str]:
    """Setup a conversion dictionary to convert codons to amino acids

    Returns:
        dict: Dictionary mapping codons of form 'xyz' to amino acids of form 'X'
    """
    # Defined bases
    bases = ["t", "c", "a", "g", "x", "z", "o"]
    # Defined amino acids in correct order
    aminoacids = (
        "FFLLXZOSSSSXZOYY!!XZOCC!WXZOXXXXXXXZZZZXZOOOOOXOOLLLLXZOPPPPXZOHHQQXZORRRRXZOXXX"
        "XXXXZZZZXZOOOOOXOOIIIMXZOTTTTXZONNKKXZOSSRRXZOXXXXXXXZZZZXZOOOOOXOOVVVVXZOAAAAXZ"
        "ODDEEXZOGGGGXZOXXXXXXXZZZZXZOOOOOXOOXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        "XXXXXZZZZXZOZZZZXZOZZZZXZOZZZZXZOXXXXXXXZZZZXZOOOOOXOOOOOOXOOOOOOXOOOOOOXOOOOOOX"
        "OOXXXXXXXOOOOXOOOOOOXOO"
    )
    all_codons = [a + b + c for a in bases for b in bases for c in bases]
    return dict(zip(all_codons, aminoacids))
