import pytest, gumpy, copy

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

    #TODO: apply_variant_file()

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
