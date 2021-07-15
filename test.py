# import cprofilev
import gumpy
import time

start = time.time()
# g1 = gumpy.Genome("config/NC_000962.3.gbk")
g1 = gumpy.Genome("config/NC_004148.2.gbk")
print("Done new: ", time.time() - start)
print(g1.contains_gene("F"))
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