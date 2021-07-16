import gumpy
import time

start = time.time()
# g1 = gumpy.Genome("config/NC_000962.3.gbk")
g1 = gumpy.Genome("config/NC_004148.2.gbk")
print("Done new: ", time.time() - start)
vcf = gumpy.VariantFile("tests/test-cases/01.vcf")
g2 = g1.apply_variant_file(vcf)
# print(g1 == g2)
# print(g1.nucleotide_sequence)
# print(g2.nucleotide_sequence)

start = time.time()
g1.save("test.json") #Dumps TB in ~13s
print("Dumped", time.time() - start)
start = time.time()
g3 = gumpy.Genome.load("test.json") #Loads TB in ~9s
print("Loaded", time.time() - start)
# print([(attr) for attr in vars(g3)])
# print([(attr) for attr in vars(g3.genes["yajC"])])
# print(g3.nucleotide_sequence)
# g1.genes["yajC"].name = "yajC - g1"
# g3.genes["yajC"].name = "yajc - g3"
# yajC = g3.genes["yajC"]
# print(yajC.nucleotide_sequence)
# print(yajC.nucleotide_number)
# print(yajC.is_cds)
# print(yajC.index)
print("Success: ", g1 == g3)
# start = time.time()
# # g2 = gumpy.Genome2("config/NC_000962.3.gbk")
# g2 = gumpy.Genome2("config/NC_004148.2.gbk")
# print("Done old: ", time.time() - start) #Longer file takes roughly 510 seconds...

# print()
# print()
# print("Genomes equal: ", g1 == g2)
# print("Genes equal: ", g1.genes == g2.genes)
# print(1, g1)
# # print(1, g1.genes.keys())
# print(2, g2)
# print(2, g2.genes.keys())

# for item in dir(g1):
#     if item[0] == "_":
#         continue
#     try:
#         print("ITEM: ", item)
#         print("New: ", getattr(g1, item))
#         print()
#         print("Old: ", getattr(g2, item))
#         print()
#         print("Same: ", getattr(g1, item) == getattr(g2, item))
#         print("************")
#         print()
#     except:
#         continue

# for (key1, key2) in zip(sorted(g1.genes.keys()), sorted(g2.genes.keys())):
#     print(key1, key2)
#     print(g1.genes[key1], g2.genes[key2])
#     print(g1.genes[key1] == g2.genes[key2])
#     print()
    
# print(dir(g1))


# g2 = gumpy.Genome("config/NC_004148.2.gbk")

# print(g1.contains_gene("N"))