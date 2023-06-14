[![Tests](https://github.com/oxfordmmm/gumpy/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/gumpy/actions/workflows/tests.yaml)
[![codecov](https://codecov.io/gh/oxfordmmm/gumpy/branch/master/graph/badge.svg)](https://codecov.io/gh/oxfordmmm/gumpy) 
[![Docs](https://github.com/oxfordmmm/gumpy/actions/workflows/docs.yaml/badge.svg)](https://oxfordmmm.github.io/gumpy/)
[![PyPI version](https://badge.fury.io/py/gumpy.svg)](https://badge.fury.io/py/gumpy)

# gumpy
Genetics with Numpy

## Installation
```
git clone https://github.com/oxfordmmm/gumpy
cd gumpy
pip install .
```
## Documentation
https://oxfordmmm.github.io/gumpy/

## Testing
A suite of tests can be run from a terminal:
```
python -m pytest --cov=gumpy -vv
```

## Usage
### Parse a genbank file
Genome objects can be created by passing a filename of a genbank file
```
from gumpy import Genome

g = Genome("filename.gbk")
```

### Parse a VCF file
VCFFile objects can be created by passing a filename of a vcf file
```
from gumpy import VCFFile

vcf = VCFFile("filename.vcf")
```

### Apply a VCF file to a reference genome
The mutations defined in a vcf file can be applied to a reference genome to produce a new Genome object containing the changes detailed in the vcf.

If a contig is set within the vcf, the length of the contig should match the length of the genome. Otherwise, if the vcf details changes within the genome range, they will be made.
```
from gumpy import Genome, VCFFile

reference_genome = Genome("reference.gbk")
vcf = VCFFile("filename.vcf")

resultant_genome = reference_genome + vcf
```

### Genome level comparisons
There are two different methods for comparing changes. One can quickly check for changes which are caused by a given VCF file. The other can check for changes between two genome. The latter is therefore suited best for comparisons in which either both genomes are mutated, or the VCF file(s) are not available. The former is best suited for cases where changes caused by a VCF want to be determined, but finding gene-level differences will require rebuilding the Gene objects, which can be time consuming.

#### Compare genomes
Two genomes of the same length can be easily compared, including equality and changes between the two.
Best suited to cases where two mutated genomes are to be compared.
```
from gumpy import Genome, GenomeDifference

g1 = Genome("filename1.gbk")
g2 = Genome("filename2.gbk")

diff = g2 - g1 #Genome.difference returns a GenomeDifference object
print(diff.snp_distance) #SNP distance between the two genomes
print(diff.variants) #Array of variants (SNPs/INDELs) of the differences between g2 and g1
```

### Gene level comparisons
When a Genome object is instanciated, it is populated with Gene objects for each gene detailed in the genbank file.
These genes can also be compared.
Gene differences can be found through direct comparison of Gene objects, or systematically through the `gene_differences()` method of `GenomeDifference`.
```
from gumpy import Genome, Gene

g1 = Genome("filename1.gbk")
g2 = Genome("filename2.gbk")

#Get the Gene objects for the gene "gene1_name" from both Genomes
g1_gene1 = g1.build_gene["gene1_name"]
g2_gene1 = g2.build_gene["gene1_name"]

g1_gene1 == g2_gene1 #Equality check of the two genes
diff= g1_gene1 - g2_gene1 #Returns a GeneDifference object
diff.mutations #List of mutations in GARC describing the variation between the two genes
```

### Save and load Genome objects
Due to how long it takes to create a Genome object, it may be beneficial to save the object to disk. The reccomendation is to utilise the `pickle` module to do so, but due to the security implications of this, do so at your own risk! An example is below:
```
import pickle

import gumpy

#Load genome
g = gumpy.Genome("filename.gbk")

#Save genome
pickle.dump(g, open("filename.pkl", "wb"))

#Load genome
g2 = pickle.load(open("filename.pkl", "rb"))

g == g2 #True
```
