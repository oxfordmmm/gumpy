import pytest, gumpy, numpy


#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def test_ribosomal_slippage():
    #Open genome with PRF
    g = gumpy.Genome("config/TEST-RNA-PRF.gbk")
    #Open genome without PRF
    ref = gumpy.Genome("config/TEST-RNA.gbk")

    #A lot of genome/gene attributes should be the same
    assert numpy.all(g.nucleotide_sequence == ref.nucleotide_sequence)
    assert numpy.all(g.genes.keys() == ref.genes.keys())
    assert g.name == ref.name
    assert numpy.all(g.nucleotide_index == ref.nucleotide_index)
    assert numpy.all(g.stacked_is_reverse_complement == ref.stacked_is_reverse_complement)

    #Checking for impact on other genes
    for gene_name in  g.genes.keys():
        if gene_name != "A":
            assert g.genes[gene_name] == ref.genes[gene_name]

    #Checking for correct handling of -1 PRF
    a = g.genes["A"]
    assert numpy.all(a.nucleotide_sequence == list("aaaaaaaaaacccccccccccggggggggggttccc"))
    assert numpy.all(a.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1,32))+[-6, -5, -4] ))



