# import cprofilev
import gumpy


g1 = gumpy.Genome("config/NC_000962.3.gbk")
print("Done new")
g2 = gumpy.Genome2("config/NC_000962.3.gbk")
print("Done old")
print(g1 == g2)
print(g1.genes == g2.genes)
print(g1)
print(g1.genes.keys())
print(g2)
print(g2.genes.keys())


# g2 = gumpy.Genome("config/NC_004148.2.gbk")

# print(g1.contains_gene("N"))