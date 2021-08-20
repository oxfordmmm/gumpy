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
    assert numpy.all(a.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1,31))+[-6, -5, -4]))
    assert numpy.all(a.index == numpy.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,97,98,99]))
    assert numpy.all(a.is_cds == numpy.array([False, False, False]+[True for i in range(1,31)]+[False, False, False]))
    assert numpy.all(a.is_promoter == numpy.array([True, True, True]+[False for i in range(1,31)]+[True, True, True]))
    assert numpy.all(a.is_indel == numpy.array([False for i in range(36)]))
    assert numpy.all(a.indel_length == numpy.array([0 for i in range(36)]))
    seq = list("aaaaaaacccccccccccggggggggggtt")
    assert numpy.all(a.codons == numpy.array([''.join(seq[i*3:i*3+3]) for i in range(10)]))
    assert numpy.all(a.amino_acid_sequence == numpy.array(list("KKTPPPGGGV")))
    assert numpy.all(a.triplet_number == numpy.array([i//3 for i in range(3,33)]))

