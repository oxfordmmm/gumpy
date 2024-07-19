"""
Classes used to parse and store VCF data
"""

import copy
import pathlib
import warnings
from collections import defaultdict
from typing import Collection, Dict, Iterable, List, Tuple
from tqdm import tqdm

import numpy
import pandas
import pysam


class VCFRecord(object):
    """
    * Class for VCF records
    * Instance variables:
        * `chrom` (str): Name of the sample.
        * `pos` (int): Genome index for the change.
        * `ref` (str): Reference value for the nucleotide.
        * `alts` (tuple(str) | None): Alternative calls. Tuple can contain single values
            and indels.
        * `qual` (None): Values for quality (values other than None have yet to be
            found in testing files...).
        * `filter` (str): Whether this record has pass the filter.
        * `info` (dict): Dictionary of key->value for the info fields.
        * `values` (dict): Dictionary of key->value for the values. Usually this is the
            FORMAT field names with their corresponding values.
        * `is_filter_pass` (bool): does the filter column contain PASS?
        * `call1` (int): the index of the first call
        * `call2` (int): the index of the second call
        * `is_reference` (bool): is the call for the reference?
        * `is_null` (bool): is the call a null call?
        * `is_heterozygous` (bool): is it a is_heterozygous call i.e. call1!=call2?
        * `is_alt` (bool): or, is the call for a single specified alt
    """

    def __init__(
        self, record: pysam.libcbcf.VariantRecord, sample: str | int, min_dp: int | None
    ):
        """Constructor for the VCFRecord object.

        Parses the supplied pysam object and presents in a more Pythonic format

        Args:
            record (pysam.libcbcf.VariantRecord): The record object
            sample (str | int) : Name of the sample to consider. Used for possible cases
                                where there is more than 1 sample per record
            min_dp (int | None): Minimum depth to consider a call.
        """

        assert (
            len(record.samples.keys()) == 1
        ), "only supporting single samples per row at present!"

        assert (
            "GT" in record.samples[sample].keys()
        ), "require GT in FORMAT column to parse genotype"

        # Save some of the easier to get to attributes
        self.chrom = record.chrom
        self.contig = record.contig
        self.pos = record.pos

        if record.ref is not None:
            self.ref = record.ref.lower()
        else:
            self.ref = ""

        self.alts: tuple | None
        if record.alts is not None:
            self.alts = tuple([i.lower() for i in record.alts])
        else:
            self.alts = None

        self.qual = record.qual

        # Get the filter attribute value
        assert len(record.filter.items()) >= 0, "A record has more than 1 filter set!"
        self.filter: str | None
        if len(record.filter.items()) == 0:
            self.filter = None
            self.is_filter_pass = False
        else:
            self.filter = str(record.filter.items()[0][0])
            self.is_filter_pass = (
                True if record.filter.items()[0][0] == "PASS" else False
            )

        # Get the info field
        self.info = {}
        for key, value in record.info.items():
            self.info[key] = value

        # Get the values
        self.values = {}
        for key, item in record.samples[sample].items():
            # incorporate the logic from the old Genotype class here
            if key == "GT":
                if len(item) == 2:
                    # Ploidy 2 is expected here
                    call1, call2 = item
                    self.call1 = call1 if call1 is not None else -1
                    self.call2 = self.call1 if call2 is None else call2
                else:
                    # GVCF edge case with ploidy 1
                    call1 = item[0]
                    self.call1 = call1 if call1 is not None else -1
                    self.call2 = self.call1
                self.is_reference = (
                    True
                    if self.call1 == 0
                    and (self.call1 == self.call2 or self.call2 == -1)
                    else False
                )
                self.is_heterozygous = True if self.call1 != self.call2 else False
                self.is_null = True if set([self.call1, self.call2]) == {-1} else False
                self.is_alt = (
                    True
                    if not self.is_reference
                    and not self.is_heterozygous
                    and not self.is_null
                    else False
                )

            # Due to how pysam reads floats, there are some erroneously long dps,
            #   so round
            elif isinstance(item, float):
                item = round(item, 3)
            self.values[key] = item
        self.values["POS"] = self.pos

        if min_dp is not None:
            allelic_depth_tag = "COV" if "COV" in self.values.keys() else "AD"
            # Ensure we have a COV tag for downstream analysis
            self.values["COV"] = self.values[allelic_depth_tag]
            if self.values[allelic_depth_tag] != (None,):
                # If the depth given is below the threshold,
                #   this row is a null call's row from the GVCF
                # So treat it as such.
                if len(self.values[allelic_depth_tag]) == 1:
                    if self.values[allelic_depth_tag][0] < min_dp:
                        # Odd case where GVCF only gives the depth of the alt
                        self.is_null = True
                        self.is_heterozygous = False
                        self.is_alt = False
                        self.is_reference = False
                elif (
                    self.values[allelic_depth_tag][self.call1] < min_dp
                    or self.values[allelic_depth_tag][self.call2] < min_dp
                ):
                    self.is_null = True
                    self.is_heterozygous = False
                    self.is_alt = False
                    self.is_reference = False

        if self.is_null:
            # Override filter fails on nulls to enforce that they are detected
            self.is_filter_pass = True

    def __repr__(self) -> str:
        """Pretty print the record

        Returns:
            str: String representation of the record
        """
        s = self.chrom + "\t"
        s += str(self.pos) + "\t"
        s += self.ref + "\t"
        s += str(self.alts) + "\t"
        if self.qual is None:
            s += ".\t"
        else:
            s += str(self.qual) + "\t"
        if self.filter is None:
            s += ".\t"
        else:
            s += self.filter + "\t"
        for val in self.values.keys():
            if val == "POS":
                continue
            s += str(val) + ":"
        s = s[:-1] + "\t"
        for key in self.values.keys():
            val = self.values[key]
            if key == "POS":
                continue
            s += str(val) + ":"
        s = s[:-1]
        s += "\n"
        return s


class VCFFile(object):
    """
    * Class to instanciate a variant file (VCF)
    * Used to apply a VCF file to a genome
    * Instance variables:
        * `filename` (str): path to the VCF file
        * `vcf_version` (tuple(int)): Tuple of ints to show the VCF version of the
            file. e.g 3.2 would be (3, 2).
        * `contig_lengths` (dict): Dictionary mapping contig_name->length for all
            defined contigs.
        * `format_fields_description` (dict): Dictionary mapping
            format_name->dict(description, id, type).
        * `records` (list(VCFRecord)): List of VCFRecord objects for each record within
            the file.
        * `calls` (dict): Dict of definite calls made in the VCF file, after any
            additional filtering has been applied
        * `ignore_filter` (bool): whether to ignore the FILTER in the VCF file
        * `format_fields_min_thresholds` (dict): dictionary specifying minimum
            thresholds to be applied to fields in the FORMAT field e.g. {'GTCONF':5}
        * `variants` (numpy.array): Numpy array of the detected variants in the VCF file
        * `nucleotide_index` (numpy.array): Array of genome indices which are affected
            by the VCF
        * `ref_nucleotides` (numpy.array): Array of REF bases
        * `alt_nucleotides` (numpy.array): Array of ALT bases
        * `indel_length` (numpy.array): Array of lengths of insertions (+ve) or
            deletions (-ve) at each site
        * `is_snp` (numpy.array): Array to act as a mask for `nucleotide_index` to
            show which are SNPs
        * `is_het` (numpy.array): Array to act as a mask for `nucleotide_index` to
            show which are heterozygous calls
        * `is_null` (numpy.array): Array to act as a mask for `nucleotide_index` to
            show which are null calls
        * `is_indel` (numpy.array): Array to act as a mask for `nucleotide_index` to
            show which are indel calls
        * `snp_distance` (int): SNP distance caused by the VCF
    """

    def __init__(
        self,
        filename: str,
        ignore_filter: bool = False,
        bypass_reference_calls: bool = False,
        format_fields_min_thresholds: Dict[str, int] | None = None,
        minor_population_indices: Collection[int] = [],
        min_dp: int | None = None,
    ):
        """
        Constructor for the VCFFile object.

        Parses the VCF file using pysam.

        Args:
            filename (str) : The name of the VCF file
            ignore_filter (bool, optional): If True, ignore the FILTER column in the
                VCF file. Default is False.
            bypass_reference_calls (bool, optional): If True, skip any row in the VCF
                (and therefore do not record)  which calls reference (i.e. 0/0).
                Default is False.
            format_fields_min_thresholds (dict, optional): Dict of field name in the
                FORMAT column and a minimum threshold to apply e.g. {'DP':5}
            minor_population_indices (set, optional): set of genome indices names
                within which to look for minor populations
            min_dp (int, optional): Minimum depth to consider a call. Default is None
        """

        self.ignore_filter = ignore_filter
        assert isinstance(self.ignore_filter, bool)

        self.bypass_reference_calls = bypass_reference_calls
        assert isinstance(self.bypass_reference_calls, bool)

        self.format_fields_min_thresholds = format_fields_min_thresholds
        if self.format_fields_min_thresholds is not None:
            assert isinstance(self.format_fields_min_thresholds, dict)

        self.minor_population_indices = minor_population_indices
        # As {}/[] is a dangerous default value, convert from None to {} as req
        if self.minor_population_indices is None:
            self.minor_population_indices = set()
        else:
            # Value given, so check if it's of the right format
            # Functionally, we don't care if it's actually a set. We just care its an
            #   interable of ints
            try:
                for i in self.minor_population_indices:
                    assert isinstance(
                        i, int
                    ), "Item in minor_population_indices is not an int: " + str(i)
            except TypeError:
                # Not iterable
                assert False, "minor_population_indices given is not iterable! " + str(
                    self.minor_population_indices
                )

        assert isinstance(filename, str)
        # Use expand user path to allow use of "~"
        self.filename = str(pathlib.Path(filename).expanduser())
        assert pathlib.Path(self.filename).is_file()

        # Use pysam to parse the VCF
        # Pylint doesn't think there is a VariantFile method but it works...
        vcf = pysam.VariantFile(self.filename)

        # Get some basic metadata
        self.vcf_version = vcf.version

        # Get the contig lengths from the header
        self.contig_lengths = {}
        for name in list(vcf.header.contigs):
            self.contig_lengths[name] = vcf.header.contigs[name].length

        # Get the formats
        self.format_fields_metadata = {}
        for format_ in vcf.header.formats.keys():
            description = vcf.header.formats[format_].description
            # mypy doesn't like this as apparently this doesn't have an 'id' attr
            #   however, it works and evidently does
            id_ = vcf.header.formats[format_].id  # type: ignore
            f_type = vcf.header.formats[format_].type
            self.format_fields_metadata[format_] = {
                "description": description,
                "id": id_,
                "type": f_type,
            }

        if isinstance(self.format_fields_min_thresholds, dict):
            assert set(self.format_fields_min_thresholds.keys()).issubset(
                set(self.format_fields_metadata.keys())
            ), "field to threshold on not found in the FORMAT column of the vcf!"

        # Get the records
        self.records = []
        for record in list(vcf):
            for sample in record.samples.keys():
                self.records.append(VCFRecord(record, sample, min_dp))

        # Find calls will ensure that no calls have same position
        self.__find_calls()

        self.__get_variants()

        self.snp_distance = numpy.sum(self.is_snp)

        if len(self.minor_population_indices) > 0:
            # We are asking to find some minor variants
            # So we need to check if the COV/AD field exists as this shows
            #   minor populations
            assert (
                "COV" in self.format_fields_metadata.keys()
                or "AD" in self.format_fields_metadata.keys()
            ), (
                "'COV' and 'AD' not in VCF format fields. "
                "No minor populations can be found!"
            )
            self._find_minor_populations()
        else:
            # Give a sensible default value otherwise
            self.minor_populations: List = []

    def __repr__(self) -> str:
        """Overload the print function to write a summary of the VCF file

        Returns:
            str: String summarising the VCF file
        """
        output = (
            "VCF variant file, version "
            + "".join(str(i) + "." for i in self.vcf_version)[:-1]
            + "\n"
        )
        output += self.filename + "\n"
        output += str(len(self.records)) + " records" + "\n"
        output += (
            "FORMAT columns: "
            + ", ".join(i for i in sorted(list(self.format_fields_metadata.keys())))
            + "\n"
            + "\n"
        )
        if len(self.records) > 3:
            output += str(self.records[0])
            output += str(self.records[1])
            output += str(self.records[2])
            output += "...\n"
            output += str(self.records[-1])
        else:
            for record in self.records:
                output += str(record) + "\n"
        return output

    def _find_minor_populations(self):
        """Find the minor populations for this VCF based on the minor population
            positions given.
        Deconstructs these into actual mutations of SNP/INDEL too as this logic already
            exists here

        * Minor populations are intentionally kept very separate from other variants.
        * They are very infrequent compared to other variants so storing in similar
            structures would be resource heavy.
        * They are instead stored consistently in the form of calls:
            * [position, type, bases, depth, frs]
                * position: Genome position (or Gene position if within genes)
                * type: String denoting type of variant. One of ['ref', 'snp',
                    'ins', 'del']
                * bases: Descripton of the bases.
                    * If type in ['ref', 'snp'], of the format (ref base, alt base)
                    * If type in ['ins', 'del'], string of the bases inserted or deleted
                * depth: number of reads supporting this call
                * frs: fractional read support (i.e depth/total coverage at this point)
        """
        self.minor_populations = []
        simple_calls = []
        for idx, type_ in self.calls.keys():
            for item in self.calls[(idx, type_)]:
                # Get the simple format of this call for comparison
                if isinstance(item["call"], tuple):
                    # Indels
                    t = item["call"][0]
                    bases = item["call"][1]
                    pos = item["pos"]
                else:
                    # Snps
                    if item["call"] == "x":
                        # Null calls shouldn't have minor populations
                        continue
                    t = "snp"
                    bases = (
                        item["ref"],
                        item["call"],
                    )
                    pos = item["pos"]
                simple_calls.append((idx, pos, t, bases))
        seen = set()

        for idx, type_ in tqdm(self.calls.keys()):
            for item in self.calls[(idx, type_)]:
                # Check if we've delt with this vcf already
                if str(item["original_vcf_row"]) in seen:
                    continue
                else:
                    seen.add(str(item["original_vcf_row"]))

                # Pull out depth tag from the specific row's format fields
                # as the file metadata isn't a guarantee of the actual
                # fields of this row
                allelic_depth_tag = (
                    "COV" if item["original_vcf_row"].get("COV", None) else "AD"
                )

                # Checking for het calls
                if item["call"] == "z":
                    if 0 not in item["original_vcf_row"]["GT"]:
                        # Het call without a wildtype call, so warn about
                        # behaviour
                        warnings.warn(
                            f"Minor population detected at position {idx}, which "
                            f"doesn't include a wildtype call. Call: "
                            f"{item['original_vcf_row']['GT']}. Note that there "
                            "may be multiple mutations given at this index",
                            UserWarning,
                        )

                # Reference base(s)
                ref = item["original_vcf_row"]["REF"]

                if item["original_vcf_row"]["ALTS"] is None:
                    # Case arrises from gvcf ref calls not giving any alts
                    calls = [item["original_vcf_row"]["REF"]]
                else:
                    # Get all of the calls
                    calls = [item["original_vcf_row"]["REF"]] + list(
                        item["original_vcf_row"]["ALTS"]
                    )

                # Break down the calls as appropriate
                simple = [self._simplify_call(ref, alt) for alt in calls]

                # Map each call to the corresponding read depth
                dps = list(item["original_vcf_row"][allelic_depth_tag])

                # GVCF null calls get None for depth, so catch (and skip) this
                if dps == [None]:
                    continue
                else:
                    total_depth = sum(dps)

                # idx here refers to the position of this call, NOT this vcf row,
                # so adjust to avoid shifting when building minor calls
                idx = idx - item["pos"]
                for calls, depth in zip(simple, dps):
                    # As we can have >1 call per simple, iter
                    for call in calls:
                        # Check that this call isn't one of the actual calls
                        if (idx + int(call[0]), *call) in simple_calls:
                            # Is an actual call so we skip
                            continue
                        # Check if there are >=2 reads to support this call
                        if depth >= 2:
                            # These are minor calls!!
                            pos = idx + int(call[0])
                            if pos not in self.minor_population_indices:
                                # We don't actually care though
                                # This has to be done here as simplifying calls can move
                                #   the position
                                continue
                            if call[1] == "snp" and call[2][0] == call[2][1]:
                                # Ref calls aren't interesting
                                continue
                            # Only tracking absolute number of reads
                            self.minor_populations.append(
                                (
                                    pos,
                                    call[1],
                                    call[2],
                                    int(depth),
                                    round(depth / total_depth, 3),
                                    item["original_vcf_row"],
                                )
                            )

    def __find_calls(self):
        """
        Private method to find changes within the genome based on the variant file.

        Creates calls dict used elsewhere.
        """

        self.calls = defaultdict(list)

        for record in self.records:
            # VCF files are 1 indexed but keep for now
            index = copy.deepcopy(record.pos)

            # if we've asked, bypass (for speed) if this is a ref call
            if self.bypass_reference_calls and record.is_reference:
                continue

            # bypass filter fails , unless we have asked to ignore filter calls
            if (
                index not in self.minor_population_indices
                and not self.ignore_filter
                and not record.is_filter_pass
            ):
                continue

            # only proceed if a dictionary has been passed (otherwise defaults to None)
            if isinstance(self.format_fields_min_thresholds, dict):
                # ok to just do since we've already check in the constructor that these
                #   fields exist in the VCF
                proceed = all(
                    record.values[i] >= self.format_fields_min_thresholds[i]
                    for i in self.format_fields_min_thresholds
                )
                if not proceed:
                    continue

            if (
                len(self.minor_population_indices) > 0
                and
                # We're looking for minor populations
                not self.ignore_filter
                and
                # We don't want to ignore the filter though
                # so to avoid picking up erroneous actual calls,
                # set this to be a ref call if it is filter fail
                not record.is_filter_pass
            ):
                # We only want to allow these through if the filter fail is MIN_FRS
                if record.filter == "MIN_FRS":
                    # Allow MIN_FRS
                    variant = record.ref
                    variant_type = "ref"
                else:
                    # Discard any other filter fails
                    continue
            elif record.is_heterozygous:
                variant = "z" * len(record.ref)
                variant_type = "het"
            elif record.is_alt:
                variant = record.alts[record.call1 - 1]
                variant_type = "snp"
            elif record.is_null:
                variant = "x" * len(record.ref)
                variant_type = "null"
            elif record.is_reference:
                variant = record.ref
                variant_type = "ref"

            # if the REF, ALT pair are the same length, check if we can decompose
            #   into SNPs
            if len(record.ref) == len(variant):
                for counter, (before, after) in enumerate(zip(record.ref, variant)):
                    metadata = {}
                    metadata["call"] = after
                    metadata["ref"] = before
                    metadata["pos"] = counter
                    vcf_info = {}
                    vcf_info = copy.deepcopy(record.values)
                    vcf_info["REF"] = record.ref
                    vcf_info["ALTS"] = record.alts
                    metadata["original_vcf_row"] = vcf_info
                    if self.calls.get((index + counter, variant_type)) is not None:
                        warnings.warn(
                            UserWarning(
                                f"Multiple calls at position {index}"
                                + f" with type {variant_type} in VCF file"
                            )
                        )
                    self.calls[(index + counter, variant_type)].append(metadata)

            # otherwise the REF, ALT pair are different lengths
            else:
                # try and simplify into a single indel and multiple SNPs
                mutations = self._simplify_call(record.ref, variant)

                for p, type_, bases in mutations:
                    # p = max(p - 1, 0)
                    if type_ in ["ins", "del"]:
                        metadata = {}
                        metadata["call"] = (type_, bases)
                        metadata["ref"] = record.ref[0]
                        metadata["pos"] = p
                        vcf_info = {}
                        vcf_info = copy.deepcopy(record.values)
                        vcf_info["REF"] = record.ref
                        vcf_info["ALTS"] = record.alts
                        metadata["original_vcf_row"] = vcf_info
                        if self.calls.get((index + p, "indel")) is not None:
                            warnings.warn(
                                UserWarning(
                                    f"Multiple calls at position {index}"
                                    + " with type indel in VCF file"
                                )
                            )
                        self.calls[(index + p, "indel")].append(metadata)

                    else:
                        metadata = {}
                        metadata["call"] = bases[1]
                        metadata["ref"] = bases[0]
                        metadata["pos"] = p
                        vcf_info = {}
                        vcf_info = copy.deepcopy(record.values)
                        vcf_info["REF"] = record.ref
                        vcf_info["ALTS"] = record.alts
                        metadata["original_vcf_row"] = vcf_info
                        if self.calls.get((index + p, variant_type)) is not None:
                            warnings.warn(
                                UserWarning(
                                    f"Multiple calls at position {index}"
                                    + f" with type {variant_type} in VCF file"
                                )
                            )
                        self.calls[(index + p, variant_type)].append(metadata)

    def _simplify_call(self, ref: str, alt: str) -> List[Tuple[int, str, str]]:
        """Private method to simplify a complex call into one indel and multiple SNPs.

        Based on finding the indel position at which there is the least SNPs

        Args:
            ref (str): Reference bases. Should match reference bases at this point
            alt (str): Alt bases.

        Returns:
            [(int, str, str)]: Returns a list of tuples of (pos, one of ['ins','del',
                'snp'], indel_bases or (ref, alt) for SNPs)
        """

        def snp_number(ref: Iterable | None, alt: Iterable | None) -> int | float:
            """Count the number of SNPs between 2 sequences

            Args:
                ref (str): Reference bases
                alt (str): Alt bases

            Returns:
                int: Number of SNPs between ref and alt
            """
            if ref is None or alt is None:
                return float("inf")
            snps = 0
            for a, b in zip(ref, alt):
                if a is not None and b is not None and a != b:
                    snps += 1
            return snps

        """The process for finding the positions for indels are almost identical
        as the process for finding a del can be interpreted as finding an ins with ref 
            and alt reversed.
        The approach here is to use a sliding window to find the position of the indel 
            where the number of SNPs is lowest.
        Assumes that there is only a single indel between the ref and the alt - there 
            may be cases where this does not work these will just produce large amounts 
            of SNPs... Could be adapted to check for multi-indel but this will scale
            exponentially the number of versions checked.
        """
        if len(ref) > len(alt):
            # Del
            length = len(ref) - len(alt)
            x = ref
            y = alt
            indel = "del"
        elif len(ref) < len(alt):
            # Ins
            length = len(alt) - len(ref)
            x = alt
            y = ref
            indel = "ins"
        else:
            # Just SNPs
            # Should only be used in minor populations
            mutations_: List = []
            for i, (r, a) in enumerate(zip(ref, alt)):
                mutations_.append([i, "snp", (r, a)])
            return mutations_
        start = 0
        current: List[str | None] = []
        current_snps = float("inf")
        mutations: List = []
        for i in range(len(y) + 1):
            y1 = (
                [y[a] for a in range(i)]
                + [None for a in range(length)]
                + [y[a] for a in range(i, len(y))]
            )
            if snp_number(y1, x) <= current_snps:
                current = y1
                current_snps = snp_number(y1, x)
                start = i
        seq = [x[i] for i in range(len(current)) if current[i] is None]
        # Add the indel
        if indel == "ins":
            # Ins after index, del at index so adjust ins
            start -= 1
        mutations.append((start, indel, "".join(seq)))
        # Check for SNPs and add those
        for i, (a, b) in enumerate(zip(x, current)):
            if a is not None and b is not None and a != b:
                if indel == "ins":
                    # PWF was (i-1) but corrected, prob
                    mutations.append((i, "snp", (b, a)))
                else:
                    mutations.append((i, "snp", (a, b)))
        return mutations

    def to_df(self) -> pandas.DataFrame:
        """Convert the VCFFile to a pandas DataFrame.

        Metadata is stored in the `attrs` attribute of the DataFrame which may break
            with some operations
        (but pandas does not currently have a robust method for metadata storage...)

        Returns:
            pandas.DataFrame: DataFrame containing all of the information from the
                VCF file
        """
        meta_data: Dict = {
            "vcf_version": self.vcf_version,
            "contig_lengths": self.contig_lengths,
            "formats": self.format_fields_metadata,
        }

        chroms = []
        pos = []
        refs = []
        alts = []
        qual = []
        infos = []
        filter_ = []
        values: Dict = {}
        for record in self.records:
            chroms.append(record.chrom)
            pos.append(record.pos)
            refs.append(record.ref)
            alts.append(record.alts)
            qual.append(record.qual)
            infos.append(record.info)
            filter_.append(record.filter)
            for key in record.values:
                if values.get(key) is None:
                    values[key] = [record.values[key]]
                else:
                    values[key].append(record.values[key])
        df = pandas.DataFrame(
            {
                "CHROM": chroms,
                "POS": pos,
                "REF": refs,
                "ALTS": alts,
                "QUAL": qual,
                "INFO": infos,
                "FILTER": filter_,
                **values,
            }
        )
        df.attrs = meta_data
        return df

    def __get_variants(self):
        """
        Private method to pull the variants out of the VCFFile object.

        Builds arrays of the variant calls, and their respective genome indices, as
            well as
        masks to show whether there is a snp, het, null or indel call at the
            corresponding genome index:
        i.e is_snp[genome.nucleotide_number == indices[i]] gives a bool to determine
            if a genome has a SNP call at this position
        """

        alts = []
        variants = []
        indices = []
        refs = []
        positions = []
        is_snp = []
        is_het = []
        is_null = []
        is_indel = []
        indel_length = []
        metadata = defaultdict(list)
        to_drop = []

        for index, type_ in sorted(list(self.calls.keys())):
            for idx, item in enumerate(self.calls[(index, type_)]):
                call = item["call"]
                alt = call
                ref = item["ref"]
                if ref == alt and self.bypass_reference_calls:
                    # If we have a call at a position which the alt is a reference call
                    # we should only care if we haven't tried to bypass ref calls
                    # Have to check here too i.e CCC->ATC will have a ref call at pos 3
                    to_drop.append((idx, (index, type_)))
                    continue
                indices.append(index)
                refs.append(ref)
                pos = item["pos"]
                positions.append(pos)
                # Update the masks with the appropriate types
                if type_ == "indel":
                    # Convert to ins_x or del_x rather than tuple
                    variant = str(index) + "_" + call[0] + "_" + str(call[1])
                    alt = call[1]
                    is_indel.append(True)
                    if call[0] == "ins":
                        indel_length.append(len(call[1]))
                    else:
                        indel_length.append(-1 * len(call[1]))
                    is_snp.append(False)
                    is_het.append(False)
                    is_null.append(False)
                elif type_ == "snp":
                    variant = str(index) + ref + ">" + call
                    is_indel.append(False)
                    indel_length.append(0)
                    is_snp.append(True)
                    is_het.append(False)
                    is_null.append(False)
                elif type_ == "het":
                    variant = str(index) + ref + ">" + call
                    is_indel.append(False)
                    indel_length.append(0)
                    is_snp.append(False)
                    is_het.append(True)
                    is_null.append(False)
                elif type_ == "null":
                    variant = str(index) + ref + ">" + call
                    is_indel.append(False)
                    indel_length.append(0)
                    is_snp.append(False)
                    is_het.append(False)
                    is_null.append(True)
                elif type_ == "ref":
                    variant = str(index) + ref + ">" + ref
                    is_indel.append(False)
                    indel_length.append(0)
                    is_snp.append(False)
                    is_het.append(False)
                    is_null.append(False)
                alts.append(alt)
                variants.append(variant)
                for key in item["original_vcf_row"]:
                    metadata[key].append(item["original_vcf_row"][key])

        # Remove ref calls as required
        for idx, key in to_drop:
            del self.calls[key][idx]
            if len(self.calls[key]) == 0:
                del self.calls[key]

        # Convert to numpy arrays for neat indexing
        self.alt_nucleotides = numpy.array(alts)
        self.variants = numpy.array(variants)
        self.nucleotide_index = numpy.array(indices)
        self.is_indel = numpy.array(is_indel)
        self.indel_length = numpy.array(indel_length)
        self.is_snp = numpy.array(is_snp)
        self.is_het = numpy.array(is_het)
        self.is_null = numpy.array(is_null)
        self.ref_nucleotides = numpy.array(refs)
        self.pos = numpy.array(positions)
        self.metadata = dict()
        for key in metadata:
            self.metadata[key] = numpy.array(metadata[key], dtype=object)
