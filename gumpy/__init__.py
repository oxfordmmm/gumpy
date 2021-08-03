#! /usr/bin/env python
'''Genetics with numpy.
Utilises numpy arrays internally for improved speed when parsing genomes and assigning promoters

Classes:
    Genome
    Gene
    VariantFile
    VCFRecord
    GenomeDifference
    Genotype
'''
# -*- coding: utf-8 -*-

__version__="0.1"

from .difference import GenomeDifference, VCFDifference, GeneDifference
from .gene import Gene
from .genotype import Genotype
from .variantfile import VariantFile
from .genome import Genome
