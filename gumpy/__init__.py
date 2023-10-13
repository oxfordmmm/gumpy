#! /usr/bin/env python3
"""Genetics with numpy.
Utilises numpy arrays internally for improved speed when parsing genomes and 
assigning promoters.

Provides features including:
    * Parsing genbank files to produce Genome (and Gene) objects;
    * Handle assignment of promoter regions to genes within a genome;
    * Parsing VCF files to produce VCFFile (and VCFRecord) objects;
    * Applying a VCF file to a reference genome to produce a new genome;
    * Finding in-depth differences between two genomes such as indels;
    * Finding in-depth differences caused by a VCF such as SNPs, het and null calls;
    * Finding in-depth differences between two genes such as SNPs, indels and 
        mutations in GARC;
    * Handle intricacies of genomes such as reverse complement genes and -1 PRF.


Classes:
    * Genome
    * Gene
    * VCFFile
    * VCFRecord
    * GenomeDifference
    * GeneDifference
"""
import importlib.metadata

# Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where MAJOR is not backwards
# compatible, but MINOR and MAINTAINANCE are
__version__ = importlib.metadata.version("gumpy")

from .difference import GenomeDifference, GeneDifference, FailedComparison  # noqa
from .gene import Gene  # noqa: F401
from .variantfile import VCFFile  # noqa: F401
from .genome import Genome  # noqa: F401
