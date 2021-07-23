import pytest, gumpy, copy, numpy, os

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

    #Testing saving and loading a genome
    #Ensure that the saves directory exists
    if not os.path.exists('tests/saves'):
        os.makedirs('tests/saves')
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

def check_eq(arr1, arr2, check):
    '''Recursive helper function to determine if 2 arrays are equal. 
    Checks all sub-arrays (if exist). Will work with list, tuple and numpy.array

    Args:
        arr1 (array-like): Array 1
        arr2 (array-like): Array 2
        check (bool): Boolean accumulator

    Returns:
        bool: True when the two arrays are equal
    '''    
    if type(arr1) != type(arr2) or check == False:
        return False
    for (e1, e2) in zip(arr1, arr2):
        if type(e1) in [list, tuple, type(numpy.array([]))]:
            check = check and check_eq(e1, e2, check)
        else:
            check = check and (e1 == e2)
    return check


def test_genome_difference():
    g1 = gumpy.Genome("config/TEST-DNA.gbk", is_reference=True)
    g2 = g1.apply_variant_file(gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf"))
    assert g1 != g2

    diff = g2 - g1
    #Default view
    assert diff.snp == 5
    assert numpy.all(diff.indices == numpy.array([2, 16, 28, 78, 90]))
    assert numpy.all(diff.nucleotides == numpy.array(["g", "t", 'z', 'z', 'x']))
    assert numpy.all(diff.codons == numpy.array(['aga', 'tcc', 'zgg', 'ttz', 'aax']))
    assert numpy.all(diff.amino_acids == numpy.array(['R', 'S', 'Z', 'Z', 'X']))
    assert numpy.all(diff.indel_indices == numpy.array([71]))
    assert numpy.all(diff.indels == numpy.array([['g', 'c', 'c']]))
    assert numpy.all(diff.het_indices == numpy.array([27, 77]))

    het_calls = numpy.array([
        [(1, '*'), (99, 't'), (100, 'c')],
        [(0, '*'), (48, numpy.array(['g', 't', 't'])), (20, 'g')]
    ], dtype=object)
    assert check_eq(diff.het_calls, het_calls, True)

    assert numpy.all(diff.mutations == numpy.array(sorted([
        'A@g-2a', 'A@S5P', 'A@Z9G', 'B@Z1G'
    ])))

    #Change the view and test all outputs
    diff.update_view("full")
    assert numpy.all(diff.nucleotides == numpy.array([
        ("g", "a"), ("t", "c"), ("z", "g"),
        ("z", "t"), ("x", "a")
    ]))
    assert numpy.all(diff.codons == numpy.array([(
        ('aga', 'aaa'), ('tcc', 'ccc'), ('zgg', 'ggg'), 
        ('ttz', 'ttt'), ('aax', 'aaa')
    )]))
    assert numpy.all(diff.amino_acids == numpy.array([
        ('R', 'K'), ('S', 'P'), ('Z', 'G'), 
        ('Z', 'F'), ('X', 'K')
    ]))
    assert check_eq(diff.indels, numpy.array([
        (numpy.array(['g', 'c', 'c']), None)
    ], dtype=object), True)

    het_calls = numpy.array([
        [[(1, '*'), (99, 't'), (100, 'c')], None],
        [[(0, '*'), (48, numpy.array(['g', 't', 't'])), (20, 'g')], None]
    ], dtype=object)
    assert check_eq(diff.het_calls, het_calls, True)



    #Testing the warning about inconsistent genes
    #So make a genome with a different name for the same gene
    g3 = copy.deepcopy(g2)
    g2.genes["D"] = g2.genes["C"]
    del g2.genes["C"]
    g2.genes_lookup["D"] = g2.genes_lookup["C"]
    del g2.genes_lookup["C"]
    g2.stacked_gene_name[g2.stacked_gene_name=="C"] = "D"
    g2._Genome__recreate_genes()#Recreate the genes

    with pytest.warns(UserWarning):
        diffd = g2 - g1
    
    #Testing cases when genomes are equal
    diff = g1 - g1
    assert diff is None

    #Testing cases when 2 different genomes are given. Neither are reference
    #This is basically just for testing mutations
    g4 = copy.deepcopy(g3)
    g3.nucleotide_sequence[91] = 'g'
    g3._Genome__recreate_genes()#Recreate the genes

    diff = g3 - g4
    diff2 = g4 - g3
    assert diff.find_mutations(g1) == ["C@A2G"]
    assert diff2.find_mutations(g1) == []
    diff.update_view("full")
    assert numpy.all(diff.find_mutations(g1) == numpy.array([
            ['C@A2G', None]
        ]))

    

   

