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
    a = g.build_gene("A")
    assert numpy.all(a.nucleotide_sequence == list("aaaaaaaaaacccccccccccggggggggggttccc"))
    assert numpy.all(a.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1,31))+[-6, -5, -4]))
    assert numpy.all(a.nucleotide_index == numpy.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,97,98,99]))
    assert numpy.all(a.is_cds == numpy.array([False, False, False]+[True for i in range(1,31)]+[False, False, False]))
    assert numpy.all(a.is_promoter == numpy.array([True, True, True]+[False for i in range(1,31)]+[True, True, True]))
    assert numpy.all(a.is_indel == numpy.array([False for i in range(36)]))
    assert numpy.all(a.indel_length == numpy.array([0 for i in range(36)]))
    seq = list("aaaaaaacccccccccccggggggggggtt")
    assert numpy.all(a.codons == numpy.array([''.join(seq[i*3:i*3+3]) for i in range(10)]))
    assert numpy.all(a.amino_acid_sequence == numpy.array(list("KKTPPPGGGV")))
    assert numpy.all(a.codon_number == numpy.array([i//3 for i in range(3,33)]))

def test_ribosomal_slippage2():
    #Open genome with PRF
    g = gumpy.Genome("config/TEST-RNA-PRF2.gbk")
    #Open genome without PRF
    ref = gumpy.Genome("config/TEST-RNA.gbk")

    #A lot of genome/gene attributes should be the same
    assert numpy.all(g.nucleotide_sequence == ref.nucleotide_sequence)
    assert numpy.all(g.genes.keys() == ref.genes.keys())
    assert g.name == ref.name
    assert numpy.all(g.nucleotide_index == ref.nucleotide_index)
    assert numpy.all(g.stacked_is_reverse_complement == ref.stacked_is_reverse_complement)

    #Checking for impact on other genes
    #Impact on gene A is expected as the promoter region is impacted by changes in C
    for gene_name in  g.genes.keys():
        if gene_name not in ["A", "C"]:
            assert g.genes[gene_name] == ref.genes[gene_name]

    #Checking for correct handling of -1 PRF
    c = g.build_gene("C")
    assert numpy.all(c.nucleotide_sequence == list("ggggggggggttttttttttaaaaaaaaaaccccccccc"))
    assert numpy.all(c.nucleotide_number == numpy.array([-i for i in range(1, 31)][::-1] + list(range(1, 10))))
    assert numpy.all(c.nucleotide_index == numpy.array(list(range(61, 95))+[94, 95, 96, 97, 98] ))
    assert numpy.all(c.is_cds == numpy.array([False for i in range(61, 91)]+[True for i in range(1, 10)]))
    assert numpy.all(c.is_promoter == numpy.array([True for i in range(61, 91)]+[False for i in range(1, 10)]))
    assert numpy.all(c.is_indel == numpy.array([False for i in range(61, 100)]))
    assert numpy.all(c.indel_length == numpy.array([0 for i in range(61, 100)]))
    assert numpy.all(c.codons == numpy.array(['ccc', 'ccc', 'ccc']))
    assert numpy.all(c.amino_acid_sequence == numpy.array(list("PPP")))
    assert numpy.all(c.codon_number == numpy.array([1, 1, 1, 2, 2, 2, 3, 3, 3]))

    #Checking for appropriate changes to A based on changes in promoter region
    a = g.build_gene("A")
    assert numpy.all(a.nucleotide_sequence == list("aaaaaaaaaaccccccccccggggggggggc"))
    assert numpy.all(a.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))+[-4]))
    assert numpy.all(a.nucleotide_index == list(range(1, 31))+[99])
    assert numpy.all(a.is_cds == [False, False, False]+[True for i in range(1, 28)]+[False])
    assert numpy.all(a.is_promoter == [True, True, True]+[False for i in range(1, 28)]+[True])
    assert numpy.all(a.is_indel == [False for i in range(1, 32)])
    assert numpy.all(a.indel_length == [0 for i in range(1, 32)])
    seq = list("aaaaaaaccccccccccgggggggggg")
    assert numpy.all(a.codons == [''.join(seq[i*3:i*3+3]) for i in range(9)])
    assert numpy.all(a.amino_acid_sequence == list("KKTPPPGGG"))
    assert numpy.all(a.codon_number == [i//3 for i in range(3,30)])
