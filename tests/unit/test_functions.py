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
    # assert gumpy.Genome.load("tests/saves/TEST-DNA.json") == g1
    #Compressed
    g1.save("tests/saves/TEST-DNA.json.gz", compression_level=1)
    # assert gumpy.Genome.load("tests/saves/TEST-DNA.json.gz") == g1

    #Saving the sequence
    g1.save_sequence("tests/saves/TEST-DNA-SEQ")
    #Reloading to check this is saved correctly
    with numpy.load("tests/saves/TEST-DNA-SEQ.npz") as seq:
        s = []
        for i in seq["sequence"]:
            s.append(i)
    assert numpy.all(g1.nucleotide_sequence == s)

    #FASTA save
    g1.save_fasta("tests/saves/TEST-DNA.fasta",fixed_length=True)
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
        "metadata for all genes/loci have been included"
    ]

def test_gene_functions():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    g1 = copy.deepcopy(genome.build_gene("A"))
    g2 = copy.deepcopy(genome.build_gene("A"))

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
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")
    g2 = g1+vcf
    assert g1 != g2
    assert numpy.all(g2.nucleotide_sequence ==  numpy.array(
                    list('axaaaxxxaactcgctgcccgzgzgzzzzgttttttttataaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc')
    ))
    diff=g1-g2
    assert numpy.all(diff.variants == numpy.array(['2a>x', '6a>x', '7a>x', '8a>x', '12c>t', '14c>g', '16c>t', '17c>g','22g>z', '24g>z', '26g>z', '27g>z', '28g>z', '29g>z', '39t>a','33_ins_tt', '37_del_t', '39_ins_g', '64_ins_ca', '73_ins_a']))
    assert numpy.all(diff.nucleotide_index==numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 39, 33, 37, 39, 64, 73]))
    assert numpy.all(diff.is_indel==numpy.array([False, False, False, False, False, False, False, False, False,False, False, False, False, False, False,  True,  True,  True,True,  True]))
    assert numpy.all(diff.is_snp==numpy.array([False, False, False, False,  True,  True,  True,  True, False,False, False, False, False, False,  True, False, False, False,False, False]))
    assert numpy.all(diff.indel_length==numpy.array([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -1, 1,  2,  1]))
    assert diff.snp_distance==5


    #Check for gene level changes
    gene_changes = []
    nucleotide_changes = []
    index_changes = []
    for gene_name in g2.genes.keys():
        gene1=g1.build_gene(gene_name)
        gene2=g2.build_gene(gene_name)
        gene_changes.append(gene2!=gene1)
        nucleotide_changes.append(numpy.any(gene1.nucleotide_sequence != gene2.nucleotide_sequence))
        index_changes.append(numpy.all(gene1.nucleotide_index == gene2.nucleotide_index))
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
    g2 = g1+gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")
    assert g1 != g2

    diff = g1-g2
    #Default view
    assert diff.snp_distance == 5
    assert numpy.all(diff.nucleotide_index == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 39, 33, 37, 39, 64, 73]))
    assert numpy.all(diff.nucleotides == numpy.array(['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z', 'a'], dtype=object))

    #Change the view and test all outputs
    diff.update_view("full",'genome')
    assert numpy.all(diff.nucleotides == numpy.array([['a', 'x'],
       ['a', 'x'],
       ['a', 'x'],
       ['a', 'x'],
       ['c', 't'],
       ['c', 'g'],
       ['c', 't'],
       ['c', 'g'],
       ['g', 'z'],
       ['g', 'z'],
       ['g', 'z'],
       ['g', 'z'],
       ['g', 'z'],
       ['g', 'z'],
       ['t', 'a']]))


def test_vcf_genetic_variation():
    #Testing the VCFFile objects' difference()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")

    # checking the variant masks
    assert numpy.all(vcf.variants == numpy.array(['2a>x', '6a>x', '7a>x', '8a>x', '12c>t', '14c>g', '16c>t', '17c>g',
       '22g>z', '24g>z', '26g>z', '27g>z', '28g>z', '29g>z', '33_ins_tt',
       '37_del_t', '39_ins_g', '39t>a', '64_ins_ca', '73_ins_a']))
    assert numpy.all(vcf.nucleotide_index == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 33, 37, 39, 39, 64, 73]))
    assert numpy.all(vcf.is_snp == [False, False, False, False,  True,  True,  True,  True, False, False, False, False, False, False, False, False, False,  True, False, False])
    assert numpy.all(vcf.is_het == [False, False, False, False, False, False, False, False, True, True, True, True, True, True, False, False, False, False, False, False])
    assert numpy.all(vcf.is_indel == [False, False, False, False, False, False, False, False, False, False, False, False, False, False,  True,  True,  True, False, True,  True])
    assert numpy.all(vcf.is_null == [True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False])

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
        'REF': numpy.array(['a', 'aaa', 'aaa', 'aaa', 'c', 'c', 'ccc', 'ccc', 'g', 'g', 'gg', 'gg', 'gg', 'gg', 't', 'tt', 'tt', 'tt', 'gg', 't'], dtype=object),

        'ALTS': numpy.array([('g',), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'), ('ggt', 'gta', 'ata'),
                ('t',), ('t', 'g'), ('tgc', 'gtg'), ('tgc', 'gtg'), ('t', 'c', 'a'),
                ('t', 'c', 'a'), ('aa', 'ct', 'at'), ('aa', 'ct', 'at'), ('aa', 't', 'a'),
                ('aa', 't', 'a'), ('ttt',), ('t',), ('agt',), ('agt',), ('cagg',), ('ta', 'at')],dtype='object')
    }

    for i in full_metadata.keys():
        assert check_eq(vcf.metadata[i], full_metadata[i], True)

    assert vcf.snp_distance == 5

    #Checking gene difference objects
    g2=g1+vcf

    for gene_name in g1.genes:

        gene1=g1.build_gene(gene_name)
        gene2=g2.build_gene(gene_name)

        g_diff = gene1-gene2

        assert isinstance(g_diff, gumpy.difference.GeneDifference)

        if gene_name=='A':
            assert numpy.all(g_diff.nucleotides==numpy.array(['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z'], dtype=object))
            assert numpy.all(g_diff.mutations==numpy.array(['K1X', 'K2X', 'T3T', 'P4R', 'P5C', 'G7Z', 'G8Z', 'G9Z', 'a-2x'], dtype='<U4'))
            assert numpy.all(g_diff.ref_nucleotides==numpy.array(['aaa', 'aaa', 'acc', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg', 'a'], dtype='<U3'))

    #Testing an edge case with 2 different indels at the same position
    # this currently fails because the VCF modifies the genome and then the REF bases in the VCF are incorrect
    # g2 = g1+gumpy.VCFFile("tests/test-cases/TEST-DNA-2.vcf")
    # diff = vcf.difference(g2)
    # assert check_eq(diff.indels, {33: 'ins_2', 37: 'del_1', 40: "ins_1", 64: 'ins_2', 73: 'ins_1'}, True)

def test_gene_difference():
    #Test the Gene.difference() method and GeneDifference() objects
    genome1 = gumpy.Genome("config/TEST-DNA.gbk")
    genome2 = genome1+gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")
    g1 = genome1.build_gene("A")
    g2 = genome2.build_gene("A")
    diff = g1-g2

    assert isinstance(diff, gumpy.GeneDifference)
    assert numpy.all(diff.nucleotides == ['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z'])
    assert numpy.all(diff.mutations == ['K1X', 'K2X', 'T3T', 'P4R', 'P5C', 'G7Z', 'G8Z', 'G9Z', 'a-2x'])
    assert numpy.all(diff.ref_nucleotides == ['aaa', 'aaa', 'acc', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg', 'a'])
    assert numpy.all(diff.amino_acid_number == [1, 2, 3, 4, 5, 7, 8, 9, None])

def test_valid_variant():
    #Test if the Gene.valid_variant() works
    genome = gumpy.Genome("config/TEST-DNA.gbk")
    gene = genome.build_gene("A")

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
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")

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
        "REF": ['a', 'a', 'aaa', 'c', 'c', 'ccc', 'g', 'g', 'gg', 'gg', 't', 'tt', 'tt', 'gg', 'gg', 't'],
        "ALTS": [('g',), ('g', 't'), ('ggt', 'gta', 'ata'), ('t',), ('t', 'g'), ('tgc', 'gtg'), ('t', 'c', 'a'), ('t', 'c', 'a'),  ('aa', 'ct', 'at'), ('aa', 't', 'a'), ('ttt',), ('t',), ('agt',), ('cagg',), ('gg',), ('ta', 'at')],
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

def test_simplify_calls():

    #Not a static function, so a VCFFile must be instanciated
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")

    assert sorted(vcf._simplify_call("ga", "c")) == sorted([(0, "snp", ("g", "c")), (1, "del", "a")])
    assert sorted(vcf._simplify_call("c", "ga")) == sorted([(0, "snp", ("c", "g")), (0, "ins", "a")])

    assert vcf._simplify_call("aaa", "a") == [(1, "del", "aa")]
    assert vcf._simplify_call("a", "aaa") == [(0, "ins", "aa")]

    assert sorted(vcf._simplify_call("acgtaa", "cgt")) == sorted(
                                                        [
                                                            (0, "snp", ("a", "c")), (1, "snp", ("c", "g")),
                                                            (2, "snp", ("g", "t")), (3, "del", "taa")
                                                        ])
    assert sorted(vcf._simplify_call("cgt", "acgtaa")) == sorted(
                                                        [
                                                            (0, "snp", ("c", "a")), (1, "snp", ("g", "c")),
                                                            (2, "snp", ("t", "g")), (2, "ins", "taa")
                                                        ])

    assert sorted(vcf._simplify_call("a", "gga")) == [(-1, "ins", "gg")]
    assert sorted(vcf._simplify_call("gga",'a')) == [(0, "del", "gg")]
