import copy
import gzip
import os

import gumpy
import numpy
import pandas
import pytest

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

    assert g1.build_genome_string(fixed_length=True,)=='aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc'
    assert g2.build_genome_string(fixed_length=True,)=='aaaaataaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc'

    g2.is_indel[7]=True
    g2.indel_length[7]=2
    g2.indel_nucleotides[7]='cc'
    assert g2.build_genome_string(fixed_length=False,)=='aaaaataccaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc'
    assert g2.build_genome_string(fixed_length=True, )=='aaaaataaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc'
    
    g2.is_indel[10]=True
    g2.indel_length[10]=-4
    # g2.indel_nucleotides[7]='cccc'
    assert g2.build_genome_string(fixed_length=True, )=='aaaaataaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc'

    #Ensure that adding verbose arg to constructor doesn't change object's value
    #It should add values to the timings dict, but this is unimportant to the values so is not checked in the __eq__
    g3 = gumpy.Genome("config/TEST-DNA.gbk", verbose=True)
    assert g1 == g3

    #Ensure that the saves directory exists
    if not os.path.exists('tests/saves'):
        os.makedirs('tests/saves')

    #Saving the sequence
    g1.save_sequence("tests/saves/TEST-DNA-SEQ")
    #Reloading to check this is saved correctly
    with numpy.load("tests/saves/TEST-DNA-SEQ.npz") as seq:
        s = []
        for i in seq["sequence"]:
            s.append(i)
    assert numpy.all(g1.nucleotide_sequence == s)

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",fixed_length=True)

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",fixed_length='yes')

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",compression='yes')

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",chars_per_line=4.0)

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",chars_per_line=-10)

    with pytest.raises(Exception) as e_info:
        g1.save_fasta("tests/saves2/TEST-DNA.fasta",nucleotides_uppercase='no')

    #FASTA save
    g1.save_fasta("tests/saves/TEST-DNA.fasta",fixed_length=True)
    #Reload FASTA
    with open("tests/saves/TEST-DNA.fasta") as f:
        data = [line.replace("\n", "") for line in f]
        header = data[0]
        data = ''.join(data[1::])
    assert header == ">TEST_DNA|TEST_DNA.1|TEST_DNA, complete genome"
    assert data == "aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc".upper()

    #Try again with compression
    g1.save_fasta("tests/saves/TEST-DNA.fasta", fixed_length=True, compression=True)
    #Reload FASTA
    with gzip.open("tests/saves/TEST-DNA.fasta.gz", 'rt') as f:
        data = [line.replace("\n", "") for line in f]
        header = data[0]
        data = ''.join(data[1::])
    assert header == ">TEST_DNA|TEST_DNA.1|TEST_DNA, complete genome"
    assert data == "aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc".upper()

    #Again with other args
    g2.save_fasta("tests/saves/TEST-DNA.fasta", fixed_length=True, nucleotides_uppercase=False, description="test description", nucleotide_index_range=(5,10))
    #Reload FASTA
    with open("tests/saves/TEST-DNA.fasta") as f:
        data = [line.replace("\n", "") for line in f]
        header = data[0]
        data = ''.join(data[1::])
    assert header == ">test description"
    assert data == "aaaaataaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc"[4:9]

    g2.save_fasta("tests/saves/TEST-DNA.fasta", fixed_length=False, nucleotides_uppercase=False, description="test description", nucleotide_index_range=(5,12))
    #Reload FASTA
    with open("tests/saves/TEST-DNA.fasta") as f:
        data = [line.replace("\n", "") for line in f]
        header = data[0]
        data = ''.join(data[1::])
    assert header == ">test description"
    assert data == "ataccaaa"


    with pytest.raises(Exception):
        g1.save_fasta("tests/save/TEST-DNA.fasta", fixed_length=True, overwrite_existing=False)


    #Len
    assert g1.length == len(g1)

    #contains_gene()
    assert g1.contains_gene("A") == True
    assert g1.contains_gene("Not_A_Gene") == False
    try:
        g1.contains_gene(None)
    except AssertionError as e:
        assert str(e) == "Gene name must be string. Gene name provided was of type: <class 'NoneType'>"
    try:
        g1.contains_gene(g1)
    except AssertionError as e:
        assert str(e) == "Gene name must be string. Gene name provided was of type: <class 'gumpy.genome.Genome'>"

    with pytest.raises(Exception) as e_info:
        g1.at_index('hello')

    with pytest.raises(Exception) as e_info:
        g1.at_index(-100)

    with pytest.raises(Exception) as e_info:
        g1.at_index(100.00)

    with pytest.raises(Exception) as e_info:
        g1.at_index(1000000)

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
        "metadata for all genes/loci have been included"
    ]

def test_gene_functions():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    with pytest.raises(Exception) as e_info:
        genome.build_gene('D')

    with pytest.raises(Exception) as e_info:
        genome.build_gene(100)

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

def test_apply_vcf():
    #Will fail if the test_instanciate.py fails for test_instanciate_vcf()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")

    # should only be able to ONLY add a vcf file
    with pytest.raises(Exception) as e_info:
        g2 = g1+'hello'

    with pytest.raises(Exception) as e_info:
        g2 = g1+"tests/test-cases/TEST-DNA.vcf"

    # check that giving a VCF that includes changes outside the genome fails
    with pytest.raises(Exception) as e_info:
        g2 = g1+gumpy.VCFFile("tests/test-cases/TEST-DNA-LONG.vcf")

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

    #Checking for gene level mutations
    b1 = g1.build_gene("B")
    b2 = g2.build_gene("B")

    diff = b1 - b2
    expected = ['G1Z', 'F4L', '6_ins_tt', '10_del_t', '12_ins_g']
    assert diff.mutations.tolist() == expected

    diff = b2 - b1
    expected = ['Z1G', 'L4F', '6_del_tt', '10_ins_t', '12_del_g']
    assert diff.mutations.tolist() == expected

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

    with pytest.raises(Exception) as e_info:
        g2 = g1-'hello'

    with pytest.raises(Exception) as e_info:
        g2 = 'hello'-g1

    diff = g1-g2
    #Default view
    assert diff.snp_distance == 5
    assert numpy.all(diff.nucleotide_index == numpy.array([ 2,  6,  7,  8, 12, 14, 16, 17, 22, 24, 26, 27, 28, 29, 39, 33, 37, 39, 64, 73]))
    assert numpy.all(diff.nucleotides == numpy.array(['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z', 'a'], dtype=object))
    assert numpy.all(diff.variants == numpy.array(['2a>x', '6a>x', '7a>x', '8a>x', '12c>t', '14c>g', '16c>t', '17c>g', '22g>z', '24g>z', '26g>z', '27g>z', '28g>z', '29g>z', '39t>a', '33_ins_tt', '37_del_t', '39_ins_g', '64_ins_ca', '73_ins_a']))
    assert numpy.all(diff.indel_length[diff.is_indel] == numpy.array([2,-1,1,2,1]))
    assert numpy.all(diff.indel_nucleotides[diff.is_indel]==numpy.array(['tt','t','g','ca','a']))

    #Change the view and test all outputs
    diff.update_view("full")
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

    # now check the other way around!
    diff2 = g2-g1
    assert numpy.all(diff2.variants == numpy.array(['2x>a', '6x>a', '7x>a', '8x>a', '12t>c', '14g>c', '16t>c', '17g>c', '22z>g', '24z>g', '26z>g', '27z>g', '28z>g', '29z>g', '39a>t', '33_del_tt', '37_ins_t', '39_del_g', '64_del_ca', '73_del_a']))
    assert numpy.all(diff2.nucleotide_index == diff.nucleotide_index)
    assert numpy.all(diff2.indel_length == -1*diff.indel_length)

def test_vcf_genetic_variation():
    #Testing the VCFFile objects' difference()
    g1 = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf", bypass_reference_calls=True)

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
        ], dtype=object),
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
            assert numpy.all(g_diff.mutations==numpy.array(['K1X', 'K2X', 'T3T', 'c9t', 'P4R', 'P5C', 'G7Z', 'G8Z', 'G9Z', 'a-2x']))
            assert numpy.all(g_diff.mutations[g_diff.is_promoter]==numpy.array(['a-2x']))
            assert numpy.all(g_diff.amino_acid_number==numpy.array([1, 2, 3, None, 4, 5, 7, 8, 9, None]))
            assert numpy.all(g_diff.nucleotide_number==numpy.array([None, None, None, 9, None, None, None, None, None, -2]))
            assert numpy.all(g_diff.gene_position==numpy.array([ 1,  2,  3, 9,  4,  5,  7,  8,  9, -2]))
            assert numpy.all(g_diff.nucleotide_index==numpy.array([None, None, None, 12, None, None, None, None, None, 2]))
            assert numpy.all(g_diff.ref_nucleotides==numpy.array(['aaa', 'aaa', 'acc', 'c', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg', 'a']))
            assert numpy.all(g_diff.alt_nucleotides==numpy.array(['aax', 'xxa', 'act', 't', 'cgc', 'tgc', 'zgz', 'gzz','zzg','x']))

        elif gene_name=='B':
            assert numpy.all(g_diff.mutations==numpy.array(['G1Z', 'F4L', '6_ins_tt', '10_del_t', '12_ins_g']))
            assert numpy.sum(g_diff.is_null)==0
            assert numpy.all(g_diff.amino_acid_number[g_diff.is_het]==numpy.array([1]))
            assert numpy.all(g_diff.amino_acid_number[g_diff.is_snp]==numpy.array([4]))
            assert numpy.all(g_diff.nucleotide_number[g_diff.is_indel]==numpy.array([6,10,12]))
            assert numpy.all(g_diff.indel_nucleotides[g_diff.is_indel] == numpy.array(['tt','t','g']))
            assert numpy.all(g_diff.indel_length[g_diff.is_indel] == numpy.array([2,-1,1]))

        elif gene_name=='C':
            assert len(g_diff.mutations)==0
    
    g3 = g1 + vcf
    del g3.genes['A']

    with pytest.warns(UserWarning):
        diff = g1 - g3
    with pytest.warns(UserWarning):
        diff = g3 - g1
    
    gene1 = g1.build_gene("A")

    gene2 = g1.build_gene("A")
    gene2.total_number_nucleotides = 76

    with pytest.raises(gumpy.FailedComparison):
        gene1 - gene2

    gene3 = g1.build_gene("A")
    gene3.name = "N"
    with pytest.warns(UserWarning):
        gene1 - gene3

    gene4 = g1.build_gene("A")
    gene4.codes_protein = False
    with pytest.raises(gumpy.FailedComparison):
        gene1 - gene4
    


            

def test_gene_difference():
    #Test the Gene.difference() method and GeneDifference() objects
    genome1 = gumpy.Genome("config/TEST-DNA.gbk")
    genome2 = genome1+gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf")
    g1 = genome1.build_gene("A")
    g2 = genome2.build_gene("A")

    with pytest.raises(Exception) as e_info:
        foo=g1+'hello'

    with pytest.raises(Exception) as e_info:
        foo=g1+2

    with pytest.raises(Exception) as e_info:
        g2 = g1+genome1.build_gene('B')

    diff = g1-g2

    assert isinstance(diff, gumpy.GeneDifference)
    assert numpy.all(diff.nucleotides == ['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z'])
    assert numpy.all(diff.mutations == ['K1X', 'K2X', 'T3T', 'c9t', 'P4R', 'P5C', 'G7Z', 'G8Z', 'G9Z', 'a-2x'])
    assert numpy.all(diff.ref_nucleotides == ['aaa', 'aaa', 'acc', 'c', 'ccc', 'ccc', 'ggg', 'ggg', 'ggg', 'a'])
    assert numpy.all(diff.amino_acid_number == [1, 2, 3, None, 4, 5, 7, 8, 9, None])

    #Try the other way around too, just to check that logic works too
    diff = g2 - g1
    assert isinstance(diff, gumpy.GeneDifference)
    assert numpy.all(diff.nucleotides == ['a', 'a', 'a', 'a', 'c', 'c', 'c', 'c', 'g', 'g', 'g', 'g', 'g', 'g'])
    assert numpy.all(diff.mutations == ['X1K', 'X2K', 'T3T', 't9c', 'R4P', 'C5P', 'Z7G', 'Z8G', 'Z9G', 'x-2a'])
    assert numpy.all(diff.ref_nucleotides == ['aax', 'xxa', 'act', 't', 'cgc', 'tgc', 'zgz', 'gzz', 'zzg', 'x'])
    assert numpy.all(diff.amino_acid_number == [1, 2, 3, None, 4, 5, 7, 8, 9, None])

    #Checking for a non-coding gene (by hacking A to be non-coding)
    genome1.genes['A']['codes_protein'] = False
    genome2.genes['A']['codes_protein'] = False

    g1 = genome1.build_gene("A")
    g2 = genome2.build_gene("A")

    diff = g1-g2
    assert isinstance(diff, gumpy.GeneDifference)
    assert numpy.all(diff.nucleotides == ['x', 'x', 'x', 'x', 't', 'g', 't', 'g', 'z', 'z', 'z', 'z', 'z', 'z'])
    assert numpy.all(sorted(diff.mutations) == sorted(['a-2x', 'a3x', 'a4x', 'a5x', 'c9t', 'c11g', 'c13t', 'c14g', 'g19z', 'g21z', 'g23z', 'g24z', 'g25z', 'g26z']))

    with pytest.warns(UserWarning) as w:
        #As this contains het calls of 1/2 etc, this should cause warnings to be raised
        genome2 = genome1+gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf", minor_population_indices=range(len(genome1)))
    assert len(w) == 3
    expected_warnings = [
        "Minor population detected at position 22, which doesn't include a wildtype call. Call: (1, 2). Note that there may be multiple mutations given at this index",
        "Minor population detected at position 26, which doesn't include a wildtype call. Call: (1, 2). Note that there may be multiple mutations given at this index",
        "Minor population detected at position 28, which doesn't include a wildtype call. Call: (1, 3). Note that there may be multiple mutations given at this index"
    ]
    for i, e in zip(w, expected_warnings):
        assert i.message.args[0] == e


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
    assert gene.valid_variant("A@K2K")
    assert gene.valid_variant("20_del_100")

    #Minority populations
    assert gene.valid_variant("A@P5C:5")
    assert gene.valid_variant("A@P5C:0.5")
    assert gene.valid_variant("A@P5C:0.05")
    assert gene.valid_variant("-2_ins_agaaat:3")
    assert gene.valid_variant("-2_ins_agaaat:375")

    #Percentage deletions
    assert gene.valid_variant("A@del_1.0")
    assert gene.valid_variant("A@del_0.0")
    assert gene.valid_variant("A@del_0.5")
    assert gene.valid_variant("A@del_0.05")
    assert gene.valid_variant("del_0.05")
    assert gene.valid_variant("del_0.05:0.02")
    assert gene.valid_variant("del_0.05:4")

    #Invalid variants
    with pytest.raises(Exception) as e_info:
        assert not gene.valid_variant("")

    with pytest.raises(Exception) as e_info:
        assert not gene.valid_variant(None)

    with pytest.raises(Exception) as e_info:
        assert not gene.valid_variant(2)

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
    assert not gene.valid_variant("a-2a")
    assert not gene.valid_variant("A@28_del_aaccggttaaccggtt")
    assert not gene.valid_variant("-2_ins_agaaat:0")
    assert not gene.valid_variant("-2_ins_agaaat:-12")
    assert not gene.valid_variant("-2_ins_agaaat:0")
    assert not gene.valid_variant("-2_ins_agaaat:0.0")
    assert not gene.valid_variant("-2_ins_agaaat:-0.2")
    assert not gene.valid_variant("-2_ins_agaaat:-")
    assert not gene.valid_variant("-2_ins_agaaat:aaa")

    with pytest.raises(AssertionError):
        gene.valid_variant(None)
    with pytest.raises(AssertionError):
        gene.valid_variant(0)
    with pytest.raises(AssertionError):
        gene.valid_variant("")
    with pytest.raises(AssertionError):
        gene.valid_variant("0")
    with pytest.raises(AssertionError):
        gene.valid_variant([1,2])
    with pytest.raises(AssertionError):
        gene.valid_variant(-10)
    with pytest.raises(AssertionError):
        gene.valid_variant(gumpy.Gene)
    with pytest.raises(AssertionError):
        gene.valid_variant("@")
        
def test_vcf_to_df():
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA.vcf", bypass_reference_calls=True)

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
        'FILTER': ['PASS', 'MIN_GT','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS','PASS'],
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
    for i in df.columns:
        assert numpy.all(data[i].equals(df[i])), 'failed on '+i
    # assert df.equals(data)

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

def test_large_deletions():
    '''Test large deletion detection
    '''
    ref = gumpy.Genome("config/TEST-DNA.gbk")
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA-4.vcf")
    sample = ref + vcf

    a = ref.build_gene("A")
    #Deletes all but first promoter
    a2 = sample.build_gene("A")
    diff = a - a2

    assert numpy.all(diff.mutations == ["-1_del_aaaaaaaaccccccccccgggggggggg", 'del_0.93'])
    assert diff.vcf_evidences == [{'GT': (1, 1), 'DP': 2, 'COV': (1, 1), 'GT_CONF': 2.05, 'REF': 'aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc', 'ALTS': ('a',)}, {'GT': (1, 1), 'DP': 2, 'COV': (1, 1), 'GT_CONF': 2.05, 'REF': 'aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc', 'ALTS': ('a',)}]

    b = ref.build_gene("B")
    #Entirely deleted
    b2 = sample.build_gene("B")
    diff = b - b2

    assert numpy.all(diff.mutations == ["del_1.0"])
    assert diff.vcf_evidences == [{'GT': (1, 1), 'DP': 2, 'COV': (1, 1), 'GT_CONF': 2.05, 'REF': 'aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc', 'ALTS': ('a',)}]

    c = ref.build_gene("C")
    #Deletes 33%, so reported as normal deletion
    c2 = sample.build_gene("C")
    diff = c - c2

    #This might look wrong, but C is revcomp
    #Bases 91, 92, 93 are deleted here
    #index:  [99 98 97 96 95 94 93 92 91]
    #number: [-3 -2 -1  1  2  3  4  5  6] --> starting at gene pos of 4

    assert numpy.all(diff.mutations == ["4_del_ggg"])
    assert diff.vcf_evidences == [{'GT': (1, 1), 'DP': 2, 'COV': (1, 1), 'GT_CONF': 2.05, 'REF': 'aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc', 'ALTS': ('a',)}]

    #Checking vcf evidence when two samples are used which affect the same base
    #Should concat the evidences
    vcf2 = gumpy.VCFFile("tests/test-cases/TEST-DNA-2.vcf")
    #Change this vcf's first entry to be base 3 rather than 2
    #This puts it in the same base as the first VCF
    vcf2.calls[(3, 'null')] = vcf2.calls[(2, 'null')]
    del vcf2.calls[(2, 'null')]

    sample2 = ref + vcf2

    diff = sample - sample2

    c1 = {key[0]: vcf.calls[key]['original_vcf_row'] for key in vcf.calls.keys()}
    c2 = {key[0]: vcf2.calls[key]['original_vcf_row'] for key in vcf2.calls.keys()}

    for idx, evidence in zip(diff.nucleotide_index, diff.vcf_evidences):
        if idx in c1.keys() and idx in c2.keys():
            #Both so concat
            assert evidence == [c1[idx], c2[idx]]
        elif idx in c1.keys():
            assert evidence == c1[idx]
        elif idx in c2.keys():
            assert evidence == c2[idx]
        else:
            assert evidence == None

            assert evidence == None
    #Repeat the entry for the start of all genes as this should be picked up in gene diff
    c1[28] = c1[3]
    c1[91] = c1[3]
    
    for gene in ref.genes.keys():
        g1 = sample.build_gene(gene)
        g2 = sample2.build_gene(gene)
        diff = g1 - g2
        for idx, evidence in zip(diff.nucleotide_index, diff.vcf_evidences):
            if idx in c1.keys() and idx in c2.keys():
                #Both so concat
                assert evidence == [c1[idx], c2[idx]]
            elif idx in c1.keys():
                assert evidence == c1[idx]
            elif idx in c2.keys():
                assert evidence == c2[idx]
            else:
                assert evidence == None


def test_misc():
    '''Misc edge case testing
    '''
    ref = gumpy.Genome("config/TEST-DNA.gbk")

    #We want to test some non-coding things here, so edit the values parsed
    ref.genes["A"]["codes_protein"] = False
    a = ref.build_gene("A")

    expected = ['A gene', '30 nucleotides', "['a' 'a' 'a']", '[-3 -2 -1]', '[]', '[]']
    assert str(a).split("\n") == expected

    #Checking gene difference
    a2 = copy.deepcopy(a)
    a2.nucleotide_sequence[2] = 't'
    diff = a - a2

    #And again for genes with no promoter
    a.is_promoter = numpy.array([False for i in a.nucleotide_sequence])
    expected = ['A gene', '30 nucleotides', 'promoter likely in adjacent gene(s)', '[]', '[]']
    assert str(a).split("\n") == expected

    #Edge cases of indels within revcomp genes

    #This VCF has 95_ins_aa and 97_del_ccc
    vcf = gumpy.VCFFile("tests/test-cases/TEST-DNA-3.vcf")
    sample = ref + vcf

    c = sample.build_gene("C")
    ref_c = ref.build_gene("C")
    diff = ref_c - c
    assert sorted(diff.mutations) == sorted(['-3_del_gg', '2_ins_tt'])
    
    #Should be the same idea from the other perspective too
    diff = c - ref_c
    assert sorted(diff.mutations) == sorted(['-3_ins_gg', '2_del_tt'])

    #Make sure gene diff can handle het calls in promoters in coding genes
    ref.genes["A"]["codes_protein"] = True
    a = ref.build_gene("A")
    a2 = copy.deepcopy(a)
    a2.nucleotide_sequence[0] = 'z'

    diff = a - a2
    assert diff.mutations == ['a-3z']

