import pytest, gumpy, copy, numpy

#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def test_genome_functions():
    '''Test all of the public functions for Genome
    '''
    #Testing equality
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    g2 = gumpy.Genome("config/TEST-DNA.gbk")
    assert g1 == g2
    g2.nucleotide_sequence[5] = "t"
    assert g1 != g2

    #TODO: Testing subtraction... Makes more sense to do this after __sub__ changes

    #Testing saving and loading a genome
    #Uncompressed
    g1.save("tests/saves/TEST-DNA.json")
    assert gumpy.Genome.load("tests/saves/TEST-DNA.json") == g1
    #Compressed
    g1.save("tests/saves/TEST-DNA.json.gz", compression_level=1)
    assert gumpy.Genome.load("tests/saves/TEST-DNA.json.gz") == g1

    #Len
    assert g1.length == len(g1)

    #contains_gene()
    assert g1.contains_gene("A") == True
    assert g1.contains_gene("Not_A_Gene") == False

    #at_index()
    try:
        g1.at_index([])
    except AssertionError as e:
        assert str(e) == "index must be an integer!"
    try:
        g1.at_index(-1)
    except AssertionError as e:
        assert str(e) == "index must be a positive integer!"
    try:
        g1.at_index(1000)
    except AssertionError as e:
        assert str(e) == "index must be less than the length of the genome!"
    
    assert g1.at_index(5) == ["A"]
    assert g1.at_index(28) == ["A", "B"]
    assert g1.at_index(65) == None

    #__repr__
    #As __repr__ returns a multiline string, convert it to a list of single line strings for ease of comparison
    r = [line for line in g1.__repr__().split("\n") if line != ""]
    assert r == [
        "TEST_DNA",
        "TEST_DNA.1",
        "TEST_DNA, complete genome",
        "99 bases",
        "aaaaaa...cccccc",
        "all genes/loci have been included"
    ]

def test_gene_functions():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    g1 = copy.deepcopy(genome.genes["A"])
    g2 = copy.deepcopy(genome.genes["A"])

    #Equality
    assert g1 == g2
    g2.name = "C"
    assert g1 != g2

    #TODO: Subtraction, changes are likely to be made to __sub__ so makes more sense to test after these changes
    
    #valid_variant() requires fixing/removing

    #__repr__
    #As __repr__ returns a multiline string, convert it to a list of single line strings for ease of comparison
    r = [line for line in g1.__repr__().split("\n") if line != ""]
    assert r == [
        "A gene",
        "30 nucleotides, codes for protein",
        "['a' 'a' 'a']",
        "[-3 -2 -1]",
        "['K' 'K' 'T' 'P' 'P' 'P' 'G' 'G' 'G']",
        "[1 2 3 4 5 6 7 8 9]"
    ]

def test_apply_vcf():
    #Will fail if the test_instanciate.py fails for test_instanciate_vcf()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")
    g2 = g1.apply_variant_file(vcf)
    assert g1 != g2
    assert numpy.all(g2.nucleotide_sequence ==  numpy.array(
                    list('agaaaaaaaaccccctccccgggggggzggttttttttttaaaaaaaaaaccccccccccggggggggggtttttttzttaaaaaaaaaxccccccccc')
    ))
    assert g2.variant_file == vcf
    assert g2.original == {
        1: 'a',
        15: 'c',
        27: 'g',
        71: 't',
        77: 't',
        89: 'a'
    }
    assert g2.indels.keys() == {71:None}.keys()
    assert numpy.all(g2.indels[71] == numpy.array(['g', 'c', 'c']))
    calls = {
        1: [(0, '*'), (68, 'g')],
        15: [(0, '*'), (68, 't')],
        27: [(1, '*'), (99, 't'), (100, 'c')],
        71: [(0, '*'), (68, numpy.array(['g', 'c', 'c']))],
        77: [(0, '*'), (48, numpy.array(['g', 't', 't'])), (20, 'g')],
        89: [(0, '*'), (68, 'x')]
    }
    assert g2.calls.keys() == calls.keys()
    for key in g2.calls.keys():
        for ((n_reads1, call1), (n_reads2, call2)) in zip(g2.calls[key], calls[key]):
            assert n_reads1 == n_reads2
            assert numpy.all(call1 == call2)
    
    #Check for gene level changes
    gene_changes = []
    nucleotide_changes = []
    index_changes = []
    for key in g2.genes.keys():
        gene_changes.append(g2.genes[key]!=g1.genes[key])
        nucleotide_changes.append(numpy.any(g2.genes[key].nucleotide_sequence != g1.genes[key].nucleotide_sequence))
        index_changes.append(numpy.all(g1.genes[key].index == g2.genes[key].index))
    assert numpy.any(gene_changes)
    assert numpy.any(nucleotide_changes)
    assert numpy.all(index_changes)