[![Build Status](https://travis-ci.com/philipwfowler/gumpy.svg?token=mdmaR7M8Ch8xBhLvrfyg&branch=master)](https://travis-ci.com/philipwfowler/gumpy) [![codecov](https://codecov.io/gh/philipwfowler/gumpy/branch/master/graph/badge.svg)](https://codecov.io/gh/philipwfowler/gumpy) [![Documentation Status](https://readthedocs.org/projects/gumpy/badge/?version=latest)](https://gumpy.readthedocs.io/en/latest/?badge=latest) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/philipwfowler/gumpy/master)

# gumpy
Genetics with Numpy

## Installation
```
git clone https://github.com/oxfordmmm/gumpy
cd gumpy
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

#### Setting a genome as a reference genome
When comparing two Genomes, it is helpful to know if a Genome is a reference genome.
This allows ease of finding mutations. To do this, use the `is_reference` flag when instanciating a Genome.
```
from gumpy import Genome

g = Genome("filename.gbk", is_reference=True)
```

### Parse a VCF file
VariantFile objects can be created by passing a filename of a vcf file
```
from gumpy import VariantFile

vcf = VariantFile("filename.vcf")
```

### Apply a VCF file to a reference genome
The mutations defined in a vcf file can be applied to a reference genome to produce a new Genome object containing the changes detailed in the vcf.

If a contig is set within the vcf, the length of the contig should match the length of the genome. Otherwise, if the vcf details changes within the genome range, they will be made.
```
from gumpy import Genome, VariantFile

reference_genome = Genome("reference.gbk")
vcf = VariantFile("filename.vcf")

resultant_genome = reference_genome.apply_variant_file(vcf)
```

### Compare genomes
Two genomes of the same length can be easily compared, including equality and changes between the two.
```
from gumpy import Genome, GenomeDifference

g1 = Genome("filename1.gbk", is_reference=True) #A reference genome
g2 = g1.apply_variant_file(VariantFile("filename2.vcf")) #A genome which has been altered by a VCF

g1 == g2 #Equality check

diff = g2 - g1 #Subtraction returns a GenomeDifference object
print(diff.snp) #SNP distance between the two genomes
print(diff.nucleotides) #Array of nucleotides in g2 which are different in g1
print(diff.indels) #Array of indels in g2 where there are indels in either g1 or g2
print(diff.het_calls) #Array of calls with coverages for every het call in both g1 and g2
print(diff.mutations) #Array of mutations within the genes in g2 compared to the reference of g1
```
### Compare VCF file
There is functionality to find the impact which a given VCF file has on a given genome.
This includes changes in codons, amino acids as well as genes
```
from gumpy import VariantFile, VCFDifference, GeneDifference, Genome
vcf = VariantFile("filename.vcf")
genome = Genome("filename.gbk", is_reference=True)

diff = vcf.difference(genome) #Returns a VCFDifference object
diff.coverage #List of the coverages of all calls

genes_diff = diff.gene_differences() #Array of GeneDifference objects
[g_diff.codons for g_diff in genes_diff] #List of the codon changes made within each gene (if the changes are within codons)
[g_diff.amino_acids for g_diff in genes_diff] #List of the amino acid changes within each gene (if coding for amino acids)
[g_diff.mutations for g_diff in genes_diff] #List of mutations within each gene in GARC
[g_diff.indels for g_diff in genes_diff] #List of the indel lengths where the indels differ between the VCF and the genome
```

### Compare genes within genomes
When a Genome object is instanciated, it is populated with Gene objects for each gene detailed in the genbank file.
These genes can also be compared.
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

