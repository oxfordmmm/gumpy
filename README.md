[![Build Status](https://travis-ci.com/philipwfowler/gumpy.svg?token=mdmaR7M8Ch8xBhLvrfyg&branch=master)](https://travis-ci.com/philipwfowler/gumpy) [![codecov](https://codecov.io/gh/philipwfowler/gumpy/branch/master/graph/badge.svg)](https://codecov.io/gh/philipwfowler/gumpy) [![Documentation Status](https://readthedocs.org/projects/gumpy/badge/?version=latest)](https://gumpy.readthedocs.io/en/latest/?badge=latest) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/philipwfowler/gumpy/master)

# gumpy
Genetics with Numpy

## Installation
```
git clone https://github.com/oxfordmmm/gumpy
cd gumpy
pip install -r requirements.txt
python setup.py build --force
pip install .
```
## Documentation
Easy access to documentation for public methods can be found using the `pydoc` module from a terminal:
```
python -m pydoc -b gumpy
```
This should open a browser window showing documentation for all loaded modules. Navigating to `gumpy (package)` should bring up available files to view documentation.

Docstrings contain documentation for almost all methods if documentation of private methods is required.

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

#### Multithreaded parsing
Due to the slow nature of some parts of instanciating a Genome object, part of it can be considerable sped up through multithreading.
Because of how multithreading works in Python, this is only beneficial when using Linux - **and so is disabled on ALL other platforms** 
To utilise mutlithreading on a Linux platform, use the `multithreaded` flag when instanciating a Genome:
```
from gumpy import Genome

g = Genome("filename.gbk", multithreaded=True)
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

resultant_genome = reference_genome.apply_variant_file(vcf)
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

diff = g2.difference(g1) #Genome.difference returns a GenomeDifference object
print(diff.snp_distance) #SNP distance between the two genomes
print(diff.nucleotides) #Array of nucleotides in g2 which are different in g1
print(diff.indels) #Array of indels in g2 where there are indels in either g1 or g2
```
#### Compare VCF file
There is functionality to find the impact which a given VCF file has on a given genome.
This includes changes in codons, amino acids as well as genes (although gene differences are computationally expensive)
```
from gumpy import VCFFile, VCFDifference, GeneDifference, Genome
vcf = VCFFile("filename.vcf")
genome = Genome("filename.gbk", is_reference=True)

diff = vcf.difference(genome) #Returns a VCFDifference object
diff.variants.get('COV') #List of the coverages of all calls
diff.snps #Dictionary mapping genome_index->snp_call

#Getting gene level differences
genes_diff = diff.gene_differences() #Array of GeneDifference objects
[g_diff.codons for g_diff in genes_diff] #List of the codon changes made within each gene (if the changes are within codons)
[g_diff.amino_acids for g_diff in genes_diff] #List of the amino acid changes within each gene (if coding for amino acids)
[g_diff.mutations for g_diff in genes_diff] #List of mutations within each gene in GARC
[g_diff.indels for g_diff in genes_diff] #List of the indel lengths where the indels differ between the VCF and the genome
```

### Gene level comparisons
When a Genome object is instanciated, it is populated with Gene objects for each gene detailed in the genbank file.
These genes can also be compared.
Gene differences can be found through direct comparison of Gene objects, or systematically through the `gene_differences()` method of both `GenomeDifference` and `VCFDifference`.
```
from gumpy import Genome, Gene

g1 = Genome("filename1.gbk")
g2 = Genome("filename2.gbk")

#Get the Gene objects for the gene "gene1_name" from both Genomes
g1_gene1 = g1.genes["gene1_name"]
g2_gene1 = g2.genes["gene1_name"]

g1_gene1 == g2_gene1 #Equality check of the two genes
g1_gene1 - g2_gene1 #Returns the indicies within the gene where the two genes differ

#Get a detailed difference between genes
diff = g1_gene1.difference(g2_gene1)
diff.mutations #List of mutations in GARC detailing the variation between the two genes
```

### Save and load Genome objects
Due to how long it takes to create a Genome object, it may be beneficial to save the object to disk.
```
from gumpy import Genome

g = Genome("filename.gbk")

g.save("output.json") #Saves the object to a JSON format

g1 = Genome.load("output.json) #Reloads the object from the JSON

g == g1 #True
```

