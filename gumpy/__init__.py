#! /usr/bin/env python3
'''Genetics with numpy.
Utilises numpy arrays internally for improved speed when parsing genomes and assigning promoters.

Provides features including:
    Parsing genbank files to produce Genome (and Gene) objects;
    Handle assignment of promoter regions to genes within a genome;
    Parsing VCF files to produce VariantFile (and VCFRecord) objects;
    Applying a VCF file to a reference genome to produce a new genome;
    Finding in-depth differences between two genomes such as indels;
    Finding in-depth differences caused by a VCF such as SNPs, het and null calls;
    Finding in-depth differences between two genes such as SNPs, indels and mutations in GARC;
    Handle intricacies of genomes such as reverse complement genes and -1 PRF.


Classes:
    Genome
    Gene
    VariantFile
    VCFRecord
    GenomeDifference
    GeneDifference
    VCFDifference
'''

__version__="0.2"

from .difference import GenomeDifference, VCFDifference, GeneDifference
from .gene import Gene
from .variantfile import VariantFile
from .genome import Genome
