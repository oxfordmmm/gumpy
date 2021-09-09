import gumpy
import time
import sys

# g = gumpy.Genome("config/NC_045512.2.gbk")
# print(g.genes_lookup)


g = gumpy.Genome("config/TEST-DNA.gbk")
# g.save("testing.json")


# x = ["" for i in list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")]
# x[28:61] = ["B" for b in range(27, 60)]
# print(len(range(27, 60)))
# x[3:31] = ["A" for a in range(4, 31)]
# x[89::] = ["C" for c in range(9)]
# print(len(range(4, 31)))
# print([j for (i, j) in enumerate(x)])
# print("".join(x))



# vcf = gumpy.VariantFile("tests/test-cases/05.vcf")
# df = vcf.to_df()
# print(df)
# print(df.attrs)
# g1 = gumpy.Genome("config/NC_004148.2.gbk", is_reference=True)
# # g2 = gumpy.Genome.load("testing.json.gz")
# # g1 = gumpy.Genome.load("test.json.gz")
# # print("Loaded...", g1==g2)
# diff = vcf.difference(g1)
# print(diff.indels)
# # print(diff.indices)
# # print(diff.coverages)
# # print(diff.het_calls)
# # print(diff.amino_acids)
# g_diff = diff.gene_differences()
# # g_diff.update_view("full")
# print([g.indel_indices.tolist() for g in g_diff])
# print([g.indels.tolist() for g in g_diff])
# print([g.mutations.tolist() for g in g_diff])
# print([g.codons.tolist() for g in g_diff])
# print([g.amino_acids.tolist() for g in g_diff])
# sys.exit()
# print([attr for attr in vars(vcf)])
# print([attr for attr in vars(vcf.records[0])])
# sys.exit()
# start = time.time()
# g1 = gumpy.Genome("config/NC_000962.3.gbk", multithreaded=True, is_reference=True)
# # # g1 = gumpy.Genome.load("test.json.gz")
# # g1 = gumpy.Genome("config/NC_004148.2.gbk", is_reference=True)
# print("Done new: ", time.time() - start)
# g2 = g1.apply_variant_file(vcf)
# diff = g1 - g2
# print(diff.snp)
# print(diff.nucleotides)
# print(diff.amino_acids)
# print(diff.indels)
# print(diff.het_calls)
# print(diff.mutations)
# # print(g1 == g2)
# # print(g1.nucleotide_sequence)
# # print(g2.nucleotide_sequence)

# start = time.time()
# g1.save("test.json.gz", compression_level=1) #Dumps TB in ~13s
# print("Dumped", time.time() - start)
# start = time.time()
# g3 = gumpy.Genome.load("test.json.gz") #Loads TB in ~9s
# print("Loaded", time.time() - start)
# print(g1 == g3)
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
# print("Success: ", g1 == g3)
# print(g1 != g2)
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