import pytest, gumpy, copy, numpy, os, pandas

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

    #Ensure that adding verbose arg to constructor doesn't change object's value
    #It should add values to the timings dict, but this is unimportant to the values so is not checked in the __eq__
    g2 = gumpy.Genome("config/TEST-DNA.gbk", verbose=True)
    assert g1 == g2

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

    #Saving the sequence
    g1.save_sequence("tests/saves/TEST-DNA-SEQ")
    #Reloading to check this is saved correctly
    with numpy.load("tests/saves/TEST-DNA-SEQ.npz") as seq:
        s = []
        for i in seq["sequence"]:
            s.append(i)
    assert numpy.all(g1.nucleotide_sequence == s)

    #FASTA save
    g1.save_fasta("tests/saves/TEST-DNA.fasta")
    #Reload FASTA
    with open("tests/saves/TEST-DNA.fasta") as f:
        data = [line.replace("\n", "") for line in f]
        header = data[0]
        data = ''.join(data[1::]).lower()
    assert header == ">TEST_DNA|TEST_DNA.1|TEST_DNA, complete genome"
    assert data == "aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc"

    #Len
    assert g1.length == len(g1)

    #contains_gene()
    assert g1.contains_gene("A") == True
    assert g1.contains_gene("Not_A_Gene") == False
    try:
        g1.contains_gene(None)
        assert False
    except AssertionError as e:
        assert str(e) == "Gene name must be string. Gene name provided was of type: <class 'NoneType'>"
    try:
        g1.contains_gene(g1)
        assert False
    except AssertionError as e:
        assert str(e) == "Gene name must be string. Gene name provided was of type: <class 'gumpy.genome.Genome'>"

    #at_index()
    try:
        g1.at_index([])
        assert False
    except AssertionError as e:
        assert str(e) == "index must be an integer!"
    try:
        g1.at_index(-1)
        assert False
    except AssertionError as e:
        assert str(e) == "index must be a positive integer!"
    try:
        g1.at_index(1000)
        assert False
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

    g2.name = "A"
    g2.nucleotide_sequence[2] = "g"
    assert g1.list_mutations_wrt(g2) == ["g-1a"]

def test_apply_vcf():
    #Will fail if the test_instanciate.py fails for test_instanciate_vcf()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")
    g2 = g1.apply_variant_file(vcf)
    assert g1 != g2
    assert numpy.all(g2.nucleotide_sequence ==  numpy.array(
                    list('axaaaxxxaactcgctgcccgzgzgzzzzgttttttttataaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc')
    ))
    assert numpy.all(g2 - g1 == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 39]))
    assert g2.variant_file == vcf
    assert g2.original == {1: 'a',
         5: 'a',
         6: 'a',
         7: 'a',
         11: 'c',
         13: 'c',
         15: 'c',
         16: 'c',
         21: 'g',
         23: 'g',
         25: 'g',
         26: 'g',
         27: 'g',
         28: 'g',
         32: 't',
         36: 't',
         38: 't',
         39: 't',
         63: 'g',
         72: 't'
    }
    assert g2.indels == {32: 'ins_2', 36: 'del_1', 39: 'ins_1', 63: 'ins_2', 72: 'ins_1'}
    calls = {
        1: [(0, '-'), (68, 'g')],
        15: [(0, '-'), (68, 't')],
        27: [(1, '-'), (99, 't'), (100, 'c')],
        71: [(0, '-'), (68, numpy.array(['g', 'c', 'c']))],
        77: [(0, '-'), (48, numpy.array(['g', 't', 't'])), (20, 'g')],
        89: [(0, '-'), (68, 'x')]
    }

    #Check for gene level changes
    gene_changes = []
    nucleotide_changes = []
    index_changes = []
    for key in g2.genes.keys():
        gene_changes.append(g2.genes[key]!=g1.genes[key])
        nucleotide_changes.append(numpy.any(g2.genes[key].nucleotide_sequence != g1.genes[key].nucleotide_sequence))
        index_changes.append(numpy.all(g1.genes[key].nucleotide_index == g2.genes[key].nucleotide_index))
    assert numpy.any(gene_changes)
    assert numpy.any(nucleotide_changes)
    assert numpy.all(index_changes)

def check_eq(arr1, arr2, check):
    '''Recursive helper function to determine if 2 arrays are equal.
    Checks all sub-arrays (if exist). Will work with list, tuple, dict and numpy.array

    Args:
        arr1 (array-like): Array 1
        arr2 (array-like): Array 2
        check (bool): Boolean accumulator

    Returns:
        bool: True when the two arrays are equal
    '''
    if not check:
        return False
    if type(arr1) == dict:
        if arr1.keys() != arr2.keys():
            return False
        else:
            for key in arr1.keys():
                e1 = arr1[key]
                e2 = arr2[key]
                if type(e1) in [list, tuple, dict, type(numpy.array([]))]:
                    check = check and check_eq(e1, e2, check)
                else:
                    check = check and (e1 == e2)
    else:
        if len(arr1) != len(arr2):
            return False
        else:
            for (e1, e2) in zip(arr1, arr2):
                if type(e1) in [list, tuple, dict, type(numpy.array([]))]:
                    check = check and check_eq(e1, e2, check)
                else:
                    check = check and (e1 == e2)
    return check


def test_genome_difference():
    g1 = gumpy.Genome("config/TEST-DNA.gbk", is_reference=True)
    g2 = g1.apply_variant_file(gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf"))
    assert g1 != g2

    diff = g2.difference(g1)
    #Default view
    assert diff.snp_distance == 15
    assert numpy.all(diff.nucleotide_index == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 39]))
    assert numpy.all(diff.nucleotides == numpy.array(['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z','z', 'a']))
    assert numpy.all(diff.indel_indices == numpy.array([33, 37, 40, 64, 73]))
    assert numpy.all(diff.indels == numpy.array(['ins_2', 'del_1', 'ins_1', 'ins_2', 'ins_1']))

    with pytest.warns(UserWarning):
        assert check_eq(diff.variants(-3), {}, True)
    with pytest.warns(UserWarning):
        assert check_eq(diff.variants(241), {}, True)
    check = True
    try:
        diff.variants(None)
        check = False #Reached if code does not assert False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw AssertationError"
    try:
        diff.variants("3")
        check = False #Reached if code does not assert False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw AssertationError"
    try:
        diff.variants([5])
        check = False #Reached if code does not assert False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw AssertationError"

    full_metadata = {
        'GT': numpy.array([(i,None) for i in [
            (None, None), (None, None), (None, None), (None, None), (1, 1), (2, 2), (1, 1), (1, 1), (1, 2), (0, 2), (1, 2),
            (1, 2), (1, 3), (1, 3), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)
        ]]),
        'DP': numpy.array([(i,None) for i in [
            2, 4, 4, 4, 50, 45, 70, 70, 202, 202, 100, 100, 100, 100, 200, 200, 200, 200, 200, 200
        ]]),
        'COV': numpy.array([(i,None) for i in [
            (1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (0, 50), (0, 2, 43), (0, 68, 8), (0, 68, 8), (1, 99, 100, 2),
            (99, 1, 100, 2), (0, 48, 50, 2), (0, 48, 50, 2), (0, 48, 2, 50), (0, 48, 2, 50), (1, 199), (1, 199), (1, 199), (1, 199), (1, 199),
            (1, 198, 1)
        ]]),
        'GT_CONF': numpy.array([(i,None) for i in [
            2.05, 2.76, 2.76, 2.76, 200.58, 155.58, 300.25, 300.25, 613.77, 613.77, 475.54, 475.54, 315.11, 315.11, 145.21, 145.21, 145.21,
            145.21, 145.21, 145.21
        ]]),
        'REF': numpy.array([(i,None) for i in [
            'a', 'aaa', 'aaa', 'aaa', 'c', 'c', 'ccc', 'ccc', 'g', 'g', 'gg', 'gg', 'gg', 'gg', 't', 'tt', 'tc', 'tc', 'gg', 't'
        ]]),
        'ALTS': numpy.array([(i,None) for i in [
            ('g',), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'), ('t',), ('t', 'g'),
            ('tgc', 'gtg'), ('tgc', 'gtg'), ('t', 'c', 'a'), ('t', 'c', 'a'), ('aa', 'ct', 'at'), ('aa', 'ct', 'at'), ('aa', 't', 'a'),
            ('aa', 't', 'a'), ('ttt',), ('t',), ('tag',), ('tag',), ('cagg',), ('ta', 'at')
        ]]),
        'call': numpy.array(
            [(i, None)
            for i in
                ['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z',
                    ('ins', 2), ('del', 1), 'a', ('ins', 1), ('ins', 2), ('ins', 1)
                ]
            ]),
        'type': numpy.array(
            [
                (i, None)
                for i in [
                    'null', 'null', 'null', 'null', 'snp', 'snp', 'snp', 'snp',
                    'het', 'het', 'het', 'het', 'het', 'het', 'indel', 'indel', 'snp', 'indel', 'indel', 'indel'
                ]
            ]),
        'pos': numpy.array([(i, None) for i in [0, 0, 1, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, -1, 0]]),
        'ref': numpy.array([(i, None) for i in ['a', 'a', 'a', 'a', 'c', 'c', 'c', 'c', 'g', 'g', 'g', 'g', 'g', 'g', 't', 't', 'c', 't', 'g', 't']])
    }
    all_indices = [2, 6, 7, 8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 33, 37, 39, 40, 64, 73]
    for (idx, genome_index) in enumerate(all_indices):
        assert check_eq(diff.variants(genome_index), {field: full_metadata[field][idx] for field in full_metadata.keys()}, True)

    #Change the view and test all outputs
    diff.update_view("full")
    assert numpy.all(diff.nucleotides == numpy.array([
       ['x', 'a'],
       ['x', 'a'],
       ['x', 'a'],
       ['x', 'a'],
       ['t', 'c'],
       ['g', 'c'],
       ['t', 'c'],
       ['g', 'c'],
       ['z', 'g'],
       ['z', 'g'],
       ['z', 'g'],
       ['z', 'g'],
       ['z', 'g'],
       ['z', 'g'],
       ['a', 't']
    ]))
    assert check_eq(diff.indels, numpy.array([['ins_2', None], ['del_1', None], ['ins_1', None], ['ins_2', None], ['ins_1', None]], dtype=object), True)


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
        diffd = g2.difference(g1)

    #Testing cases when genomes are equal
    diff = g1.difference(g1)
    assert diff is None


def test_vcf_difference():
    #Testing the VariantFile objects' difference()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")

    #Get the difference
    diff = vcf.difference(g1)
    assert isinstance(diff, gumpy.VCFDifference)
    #Checking the variant masks
    assert numpy.all(diff.variants == ['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z', 'ins_2', 'del_1', 'a', 'ins_1', 'ins_2', 'ins_1'])
    assert numpy.all(diff.nucleotide_index == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 33, 37, 39, 40, 64, 73]))
    assert numpy.all(diff.is_snp == [False, False, False, False, True, True, True, True, False, False, False, False, False, False, False, False, True, False, False, False])
    assert numpy.all(diff.is_het == [False, False, False, False, False, False, False, False, True, True, True, True, True, True, False, False, False, False, False, False])
    assert numpy.all(diff.is_indel == [False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, True, False, True, True, True])
    assert numpy.all(diff.is_null == [True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False])

    full_metadata = {
        'GT': numpy.array([
            (None, None), (None, None), (None, None), (None, None), (1, 1), (2, 2), (1, 1), (1, 1), (1, 2), (0, 2), (1, 2),
            (1, 2), (1, 3), (1, 3), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)
        ]),
        'DP': numpy.array([
            2, 4, 4, 4, 50, 45, 70, 70, 202, 202, 100, 100, 100, 100, 200, 200, 200, 200, 200, 200
        ]),
        'COV': numpy.array([
            (1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (0, 50), (0, 2, 43), (0, 68, 8), (0, 68, 8), (1, 99, 100, 2),
            (99, 1, 100, 2), (0, 48, 50, 2), (0, 48, 50, 2), (0, 48, 2, 50), (0, 48, 2, 50), (1, 199), (1, 199), (1, 199), (1, 199), (1, 199),
            (1, 198, 1)
        ]),
        'GT_CONF': numpy.array([
            2.05, 2.76, 2.76, 2.76, 200.58, 155.58, 300.25, 300.25, 613.77, 613.77, 475.54, 475.54, 315.11, 315.11, 145.21, 145.21, 145.21,
            145.21, 145.21, 145.21
        ]),
        'REF': numpy.array([
            'a', 'aaa', 'aaa', 'aaa', 'c', 'c', 'ccc', 'ccc', 'g', 'g', 'gg', 'gg', 'gg', 'gg', 't', 'tt', 'tc', 'tc', 'gg', 't'
        ]),
        'ALTS': numpy.array([
            ('g',), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'), ('t',), ('t', 'g'),
            ('tgc', 'gtg'), ('tgc', 'gtg'), ('t', 'c', 'a'), ('t', 'c', 'a'), ('aa', 'ct', 'at'), ('aa', 'ct', 'at'), ('aa', 't', 'a'),
            ('aa', 't', 'a'), ('ttt',), ('t',), ('tag',), ('tag',), ('cagg',), ('ta', 'at')
        ]),
    }
    assert check_eq(diff.metadata, full_metadata, True)

    with pytest.warns(UserWarning):
        assert check_eq(diff.variants_by_index(0), {}, True)
    with pytest.warns(UserWarning):
        assert check_eq(diff.variants_by_index(912), {}, True)
    with pytest.warns(UserWarning):
        assert check_eq(diff.variants_by_index(-10), {}, True)
    check = True
    try:
        diff.variants_by_index(None)
        check = False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw expected AssertationError"

    try:
        diff.variants_by_index([1])
        check = False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw expected AssertationError"
    assert check_eq(diff.variants_by_index(1), {}, True)
    #Check all indices
    for (idx, genome_index) in enumerate(diff.nucleotide_index.tolist()):
        assert check_eq(diff.variants_by_index(genome_index), {field: full_metadata[field][idx] for field in full_metadata.keys()}, True)

    assert diff.snp_distance == 15
    assert check_eq(diff.snps, {
        2: 'a>x', 6: 'a>x',
        7: 'a>x',  8: 'a>x',
        12: 'c>t', 14: 'c>g',
        16: 'c>t', 17: 'c>g',
        22: 'g>z', 24: 'g>z',
        26: 'g>z', 27: 'g>z',
        28: 'g>z', 29: 'g>z',
        39: 't>a'
    }, True)
    assert check_eq(diff.indels, {33: 'ins_2', 37: 'del_1', 40: 'ins_1', 64: 'ins_2', 73: 'ins_1'}, True)

    #Checking gene difference objects
    g_diff = diff.gene_differences()
    assert numpy.all([isinstance(g, gumpy.GeneDifference) for g in g_diff])
    assert check_eq([g.nucleotides for g in g_diff], [
        numpy.array(['a','a','a','a','c','c','c','c','g','g','g','g','g','g']),
        numpy.array(['g','g', 't']),
        numpy.array([])
    ], True)
    assert check_eq([g.mutations for g in g_diff], [
        numpy.array(sorted(['A@G7Z','A@G8Z', 'A@G9Z', 'A@K1X', 'A@K2X','A@P4R', 'A@P5C', 'A@3=', 'A@a-2x'])),
        numpy.array(sorted(['B@6_ins_2', 'B@10_ins_1', 'B@G1Z', 'B@13_ins_1', 'B@F4L'])),
        numpy.array([])
    ], True)
    assert check_eq([g.indel_indices for g in g_diff], [
        numpy.array([]), numpy.array([6,10, 13]), numpy.array([])
    ], True)
    assert check_eq([g.codons for g in g_diff], [
        numpy.array(['aaa', 'aaa', 'acc', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg']),
        numpy.array(["ggg", "ttt"]),
        numpy.array([])
    ], True)
    assert check_eq([g.amino_acids for g in g_diff], [
        numpy.array(['K', 'K' ,'P' ,'P' ,'G' ,'G' ,'G']),
        numpy.array(['G', 'F']),
        numpy.array([])
    ], True)

    assert check_eq(g_diff[0].nucleotide_variants(3), {
        'ALTS': (None, ('ggt', 'gta', 'ata')),
        'REF': (None, 'aaa'),
        'call': (None, 'x'),
        'type': (None, 'null')
        }, True)
    assert check_eq(g_diff[2].nucleotide_variants(1), {}, True)
    with pytest.warns(UserWarning):
        assert g_diff[0].nucleotide_variants(100) == {}
    with pytest.warns(UserWarning):
        assert g_diff[0].nucleotide_variants(0) == {}
    with pytest.warns(UserWarning):
        assert g_diff[0].nucleotide_variants(-8) == {}
    check = True
    try:
        g_diff[0].nucleotide_variants(None)
        check = False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw expected AssertationError"

    assert check_eq(g_diff[0].amino_acid_variants(1), {
        'ALTS': ([None], [('ggt', 'gta', 'ata')]),
        'REF': ([None], ['aaa']),
        'call': ([None], ['x']),
        'type': ([None], ['null'])
    }, True)
    assert check_eq(g_diff[1].amino_acid_variants(1), {
        'REF': ([None, None], ['gg', 'gg']),
        'ALTS': ([None, None], [('aa', 't', 'a'), ('aa', 't', 'a')]),
        'call': ([None, None], ['z', 'z']),
        'type': ([None, None], ['het', 'het'])
    }, True)
    assert check_eq(g_diff[2].amino_acid_variants(1), {}, True)
    assert check_eq(g_diff[1].amino_acid_variants(3), {}, True)
    with pytest.warns(UserWarning):
        assert g_diff[1].amino_acid_variants(-6) == {}
    with pytest.warns(UserWarning):
        assert g_diff[1].amino_acid_variants(0) == {}
    with pytest.warns(UserWarning):
        assert g_diff[1].amino_acid_variants(1824) == {}
    check = True
    try:
        g_diff[0].amino_acid_variants(None)
        check = False
    except AssertionError:
        pass
    finally:
        if not check:
            assert False, "Code did not throw expected AssertationError"

    #Testing an edge case with 2 different indels at the same position
    g2 = g1.apply_variant_file(gumpy.VariantFile("tests/test-cases/TEST-DNA-2.vcf"))
    diff = vcf.difference(g2)
    assert check_eq(diff.indels, {33: 'ins_2', 37: 'del_1', 40: "ins_1", 64: 'ins_2', 73: 'ins_1'}, True)

def test_gene_difference():
    #Test the Gene.difference() method and GeneDifference() objects
    genome1 = gumpy.Genome("config/TEST-DNA.gbk")
    genome2 = genome1.apply_variant_file(gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf"))
    g1 = genome1.genes["A"]
    g2 = genome2.genes["A"]
    diff = g1.difference(g2)

    assert isinstance(diff, gumpy.GeneDifference)
    assert numpy.all(diff.nucleotides == ['a', 'a', 'a', 'a', 'c', 'c', 'c', 'c', 'g', 'g', 'g', 'g', 'g', 'g'])
    assert numpy.all(diff.mutations == sorted(['A@a-2x', 'A@K1X', 'A@K2X', 'A@3=', 'A@P4R', 'A@P5C', 'A@G7Z', 'A@G8Z', 'A@G9Z']))
    assert diff.indel_indices.tolist() == []
    assert numpy.all(diff.codons == ['aaa', 'aaa', 'acc', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg'])
    assert numpy.all(diff.amino_acids == ['K', 'K' ,'P' ,'P' ,'G' ,'G' ,'G'])

def test_valid_varaint():
    #Test if the Gene.valid_variant() works
    genome = gumpy.Genome("config/TEST-DNA.gbk")
    gene = genome.genes["A"]

    #Some normal valid varaints
    assert gene.valid_variant("K2X")
    assert gene.valid_variant("K2P")
    assert gene.valid_variant("P5C")
    assert gene.valid_variant("a-2z")
    assert gene.valid_variant("a-3g")
    assert gene.valid_variant("3=")
    assert gene.valid_variant("17_ins_a")
    assert gene.valid_variant("12_ins_1")
    assert gene.valid_variant("23_ins")
    assert gene.valid_variant("-2_ins_agaaat")
    assert gene.valid_variant("12_indel")
    assert gene.valid_variant("4_del")
    assert gene.valid_variant("5_del_aa")
    assert gene.valid_variant("7_del_3")
    assert gene.valid_variant("a3c")

    assert gene.valid_variant("A@P5C")
    assert gene.valid_variant("A@a-2z")
    assert gene.valid_variant("A@8_indel")
    assert gene.valid_variant("A@9=")

    #Invalid variants
    assert not gene.valid_variant("A2X")
    assert not gene.valid_variant("A192X")
    assert not gene.valid_variant("B@K2X")
    assert not gene.valid_variant("52=")
    assert not gene.valid_variant("A@a-19c")
    assert not gene.valid_variant("aaaaaa")
    assert not gene.valid_variant("18-")
    assert not gene.valid_variant("  ")
    assert not gene.valid_variant("This is not a variant")
    assert not gene.valid_variant("A@5_del_gg")
    assert not gene.valid_variant("17_indel_2")
    assert not gene.valid_variant("B@4_ins_a")
    assert not gene.valid_variant(" @K2A")
    assert not gene.valid_variant("gene@a4t")
    assert not gene.valid_variant("A@K2K")
    assert not gene.valid_variant("a-2a")
    assert not gene.valid_variant("A@28_del_aaccggttaaccggtt")
    assert not gene.valid_variant("20_del_100")

    def assert_throws(mutation):
        check = True
        try:
            gene.valid_variant(mutation)
            check = False
        except AssertionError:
            pass
        finally:
            if not check:
                assert False, "Code did not throw expected AssertationError"

    assert_throws(None)
    assert_throws(0)
    assert_throws("")
    assert_throws("0")
    assert_throws([1,2])
    assert_throws(-10)
    assert_throws(gumpy.Gene)
    assert_throws("@")

def test_vcf_to_df():
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")

    df = vcf.to_df()
    assert df.attrs == {
        "vcf_version": (4, 2),
        "contig_lengths": {"TEST_DNA": 99},
        "formats": {
            "COV": {
                "id": 1,
                "description": "Number of reads on ref and alt alleles",
                "type": "Integer"
            },
            "GT": {
                "id": 2,
                "description": "Genotype",
                "type": "String"
            },
            "DP": {
                "id": 3,
                "description": "total kmer depth from gramtools",
                "type": "Integer"
            },
            "GT_CONF": {
                "id": 4,
                "description": "Genotype confidence. Difference in log likelihood of most likely and next most likely genotype",
                "type": "Float"
            }
        }
    }
    #Building the reference dataframe
    data = {
        "CHROM": ["TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA","TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA", "TEST_DNA"],
        "POS": [2, 4, 6, 12, 14, 16, 22, 24, 26, 28, 33, 36, 39, 65, 69, 73],
        "REF": ['a', 'a', 'aaa', 'c', 'c', 'ccc', 'g', 'g', 'gg', 'gg', 't', 'tt', 'tc', 'gg', 'gg', 't'],
        "ALTS": [('g',),
                 ('g', 't'),
                 ('ggt', 'gta', 'ata'),
                 ('t',),
                 ('t', 'g'),
                 ('tgc', 'gtg'),
                 ('t', 'c', 'a'),
                 ('t', 'c', 'a'),
                 ('aa', 'ct', 'at'),
                 ('aa', 't', 'a'),
                 ('ttt',),
                 ('t',), ('tag',),
                 ('cagg',), ('gg',), ('ta', 'at')],
        "QUAL": [None, None, None, None, None, None,None, None, None, None, None, None, None, None, None, None],
        "INFO": [{"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15},{"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {"KMER": 15}, {'KMER': 15}, {'KMER': 15}, {'KMER': 15}, {'KMER': 15}],
        "GT": [(None, None),
             (None, None),
             (None, None),
             (1, 1),
             (2, 2),
             (1, 1),
             (1, 2),
             (0, 2),
             (1, 2),
             (1, 3),
             (1, 1),
             (1, 1),
             (1, 1),
             (1, 1),
             (0, 0),
             (1, 1)],
        "DP": [2, 4, 4, 50, 45, 70, 202, 202, 100, 100, 200, 200, 200, 200, 200, 200],
        "COV": [(1, 1),
                 (1, 2, 2),
                 (1, 1, 1, 1),
                 (0, 50),
                 (0, 2, 43),
                 (0, 68, 8),
                 (1, 99, 100, 2),
                 (99, 1, 100, 2),
                 (0, 48, 50, 2),
                 (0, 48, 2, 50),
                 (1, 199),
                 (1, 199),
                 (1, 199),
                 (1, 199),
                 (1, 199),
                 (1, 198, 1)],
        "GT_CONF": [2.05,
                     3.77,
                     2.76,
                     200.58,
                     155.58,
                     300.25,
                     613.77,
                     613.77,
                     475.54,
                     315.11,
                     145.21,
                     145.21, 145.21, 145.21, 145.21, 145.21]
    }
    #Already tested the metadata so test equality against just data
    data = pandas.DataFrame(data)
    data.attrs = {}
    df.attrs = {}
    assert df.equals(data)

def test_indels():
    #Tests the VariantFile.indel()

    #Not a static function, so a VariantFile must be instanciated
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")
    assert sorted(vcf.simplify_call("ga", "c")) == sorted([(0, "snp", ("g", "c")), (1, "del", "a")])
    assert vcf.simplify_call("aaa", "a") == [(1, "del", "aa")]
    assert sorted(vcf.simplify_call("acgtaa", "cgt")) == sorted(
                                                        [
                                                            (0, "snp", ("a", "c")), (1, "snp", ("c", "g")),
                                                            (2, "snp", ("g", "t")), (3, "del", "taa")
                                                        ])
    assert sorted(vcf.simplify_call("a", "gga")) == [(-1, "ins", "gg")]
