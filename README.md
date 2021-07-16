[![Build Status](https://travis-ci.com/philipwfowler/gumpy.svg?token=mdmaR7M8Ch8xBhLvrfyg&branch=master)](https://travis-ci.com/philipwfowler/gumpy) [![codecov](https://codecov.io/gh/philipwfowler/gumpy/branch/master/graph/badge.svg)](https://codecov.io/gh/philipwfowler/gumpy) [![Documentation Status](https://readthedocs.org/projects/gumpy/badge/?version=latest)](https://gumpy.readthedocs.io/en/latest/?badge=latest) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/philipwfowler/gumpy/master)

# gumpy
Genetics with Numpy

## Usage
### Parse a genbank file
Genome objects can be created by passing a filename of a genbank file
```
from gumpy import Genome
g = Genome("filename.gbk")
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
vcf = VariantFile("filename.vcf)
resultant_genome = reference_genome.apply_variant_file(vcf)
```

### Compare genomes
Two genomes of the same length can be easily compared, including equality and changes between the two.
```
from gumpy import Genome
g1 = Genome("filename1.gbk")
g2 = Genome("filename2.gbk")
g1 == g2 #Equality check
g1 - g2 #Returns a numpy array of the indicies where the two genomes differ in nucleotide sequences
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

