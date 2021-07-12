# import cprofilev
import gumpy


g1 = gumpy.Genome("config/NC_045512.2.gbk")
print("Done new")
g2 = gumpy.Genome2("config/NC_045512.2.gbk")
print("Done old")
print()
print()
print(g1 == g2)
print(g1.genes == g2.genes)
print(1, g1)
print(1, g1.genes.keys())
print(2, g2)
print(2, g2.genes.keys())


# g2 = gumpy.Genome("config/NC_004148.2.gbk")

# print(g1.contains_gene("N"))