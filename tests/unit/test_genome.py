import pytest, numpy, copy

from pathlib import Path

from gumpy.genome import Genome

TEST_CASE_DIR = "tests/test-cases/"

reference=Genome(genbank_file="config/NC_004148.2.gbk",name="HMPV")

def test_Genome_instantiate_genbank():

    # check that the M. tuberculosis H37rV genome is the right length
    assert reference.genome_length==13335

    # check the species is stored correctly
    assert reference.organism=='Human metapneumovirus'

    # check that the sequence starts and ends as we expect
    assert reference.genome_coding_strand[0]=='a'
    assert reference.genome_coding_strand[-1]=='t'

    assert len(reference.gene_names)==9

    # try a normal gene
    mask=reference.genome_feature_name=="M2"

    sequence=reference.genome_coding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tta"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="aca"

    sequence=reference.genome_noncoding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="aat"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="tgt"

    sequence=reference.genome_sequence[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tta"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="aca"

    m2_sequence="".join(sequence)
    assert m2_sequence=='ttaattaaaaataaaataaaatttgggacaaatcataatgtctcgcaaggctccatgcaaatatgaagtgcggggcaaatgcaacagaggaagtgagtgtaagtttaaccacaattactggagttggccagatagatacttattaataagatcaaactatctattaaatcagcttttaaggaacactgatagagctgatggcctatcaataatatcaggcgcaggcagagaagacagaacgcaagattttgttctaggttccaccaatgtggttcaaggttatattgatgataaccaaagcataacaaaagctgcagcctgctacagtctacacaacataatcaagcaactacaagaagttgaagttaggcaggctagagatagcaaactatctgacagcaagcatgtggcactccataacttaatcttatcttacatggagatgagcaaaactcccgcatctttaatcaacaatctcaaaagactgccgagagaaaaactgaaaaaattagcaaagctgataattgacttatcagcaggcgctgaca'

    # try a normal gene
    mask=reference.genome_feature_name=="M2_2"

    sequence=reference.genome_coding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="atg"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="tag"

    sequence=reference.genome_noncoding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tac"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="atc"

    sequence=reference.genome_sequence[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="atg"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="tag"

    m2_2_sequence="".join(sequence)

    assert m2_2_sequence=='atgactcttcatatgccctgcaagacagtgaaagcattaatcaagtgcagtgagcatggtcctgttttcattactatagaggttgatgaaatgatatggactcaaaaagaattaaaagaagctttgtccgatgggatagtgaagtctcacaccaacatttacaattgttatttagaaaacatagaaattatatatgtcaaggcttacttaagttag'


reference2=Genome(fasta_file="config/NC_004148.2.fasta.gz",name="HMPV")

def test_Genome_instantiate_fasta():

    # check that the M. tuberculosis H37rV genome is the right length
    assert reference2.genome_length==13335

    # check the species is stored correctly
    assert reference2.organism=='Human metapneumovirus'

    # check that the sequence starts and ends as we expect
    assert reference2.genome_coding_strand[0]=='a'
    assert reference2.genome_coding_strand[-1]=='t'

def test_Genome_valid_gene_mutation_snps():

    # correct protein SNPs
    assert reference.valid_gene_mutation("F@M1N")
    assert reference.valid_gene_mutation("F@S2?")
    assert reference.valid_gene_mutation("F@S2=")
    assert reference.valid_gene_mutation("F@S2!")
    assert reference.valid_gene_mutation("M2@T76L")
    assert reference.valid_gene_mutation("M2@*?")
    assert reference.valid_gene_mutation("M2@-*?")

    # just badly formed
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("____")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("o_o_o_o_9")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("flkgjslkjg")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("o_i9k")

    # genes not present
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("lkdfjlksdjf_P1N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("rpoB@P76L")

    # incorrect reference amino acids
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@P1N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("M2@P76L")

    # bad reference amino acids
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@;1N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@B1N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@81N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("M2@J76L")

    # bad positions
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@PKN")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@P-2N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@P1000N")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@P:N")

    # bad target amino acids
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@P1O")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@PKB")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@PK;")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@PKJ")


def test_Genome_valid_gene_mutation_indels():

    # correct INDELs with good grammar
    assert reference.valid_gene_mutation("F@1_indel")
    assert reference.valid_gene_mutation("F@1_ins")
    assert reference.valid_gene_mutation("F@1_del")
    assert reference.valid_gene_mutation("F@1_ins_3")
    assert reference.valid_gene_mutation("F@1_ins_ctga")
    assert reference.valid_gene_mutation("F@1_fs")
    assert reference.valid_gene_mutation("F@1_del_3")
    assert reference.valid_gene_mutation("F@-1_indel")
    assert reference.valid_gene_mutation("F@1_del_acgt")


    # bad grammar
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@1_indl")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@1_frameshift")
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@1_ins_ggaf")

    # wrong ordering
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@indel_1")

    # incorrect gene
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F1@1_indel")

    # not in gene
    with pytest.raises(Exception):
        assert reference.valid_gene_mutation("F@2000_indel")

def test_Genome_valid_genome_variant():

    assert reference.valid_genome_variant("1a>c")
    assert reference.valid_genome_variant("2c>g")
    assert reference.valid_genome_variant("3g>a")
    assert reference.valid_genome_variant("13335t>g")

    # incorrect reference base
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("1t>c")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("2a>g")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("3t>a")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("13335c>g")

    # badly formed reference base
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("11c")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("2?>g")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("3P>a")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("13335 >g")

    # out of range index
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("0a>c")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("-1c>g")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("-2g>a")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("13336t>g")

    # badly formed index
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("1.1a>c")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("2c>fg")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("yg>a")
    with pytest.raises(Exception):
        assert reference.valid_genome_variant("tg")

def test_Genome_convert_variant_to_mutation():

    assert reference.convert_variant_to_mutation("N@a1g")=='N@M1V'
    assert reference.convert_variant_to_mutation("N@a1c")=='N@M1L'
    assert reference.convert_variant_to_mutation("N@a1t")=='N@M1L'
    assert reference.convert_variant_to_mutation("N@t2c")=='N@M1T'
    assert reference.convert_variant_to_mutation("N@t2a")=='N@M1K'
    assert reference.convert_variant_to_mutation("N@t2g")=='N@M1R'

    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@a1g")=='N@M1L'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N1@a1g")=='N@M1V'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@a1g_3")=='N@M1L'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@t1g")=='N@M1L'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@y1g")=='N@M1L'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@a1o")=='N@M1L'
    with pytest.raises(Exception):
        assert reference.convert_variant_to_mutation("N@a-1g")=='N@M1L'

def test_Genome_gbk_fasta_identical():

    assert reference.genome_length==reference2.genome_length

    assert numpy.array_equal(reference.genome_coding_strand,reference2.genome_coding_strand)

    assert numpy.array_equal(reference.genome_index,reference2.genome_index)


def test_Genome___repr__():

    assert reference.__repr__()=='NC_004148.2\nHuman metapneumovirus\nHMPV\n13335 bases\nacg...cgt'


def test_Genome___sub__():

    sample=copy.deepcopy(reference)

    sample.genome_coding_strand[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based

    indices=reference-sample

    assert indices[0]==3

    assert sample.genome_coding_strand[2]=='t'

def test_Genome_contains_gene():

    assert reference.contains_gene("F")==True
    assert reference.contains_gene("FG")==False
    assert reference.contains_gene(5)==False
    assert reference.contains_gene("")==False
    assert reference.contains_gene("rpoBC")==False
    assert reference.contains_gene("RPOB")==False

def test_Genome_at_index():

    assert reference.at_index(4686)=='F'
    assert reference.at_index(4687)=='M2'
    assert reference.at_index(5293)=='M2_2'
    assert reference.at_index(5450)=='M2_2'
    assert reference.at_index(5451)=='SH'
    assert reference.at_index(7032) is None
    assert reference.at_index(7033)=="L"
    assert reference.at_index(13150)=="L"
    assert reference.at_index(13151) is None

    # wrong gene
    assert reference.at_index(7033)!="F"
    assert reference.at_index(7033) is not None

    # bad position
    with pytest.raises(Exception):
        reference.at_index(-2)
    with pytest.raises(Exception):
        reference.at_index(0)
    with pytest.raises(Exception):
        reference.at_index(1.3)
    with pytest.raises(Exception):
        reference.at_index('gh')
    with pytest.raises(Exception):
        reference.at_index(13336)

def test_Genome_calculate_snp_distance():

    sample=copy.deepcopy(reference)

    sample.genome_coding_strand[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based
    assert sample.snp_distance(reference)==1

    # reverse the change
    sample.genome_coding_strand[2]='g'
    assert sample.snp_distance(reference)==0

    # now change two bases
    sample.genome_coding_strand[2]='t'
    sample.genome_coding_strand[3]='t'
    assert sample.snp_distance(reference)==2

sample_01=copy.deepcopy(reference)
sample_01.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF','GT_CONF_PERCENTILE'],total_coverage_threshold=5,metadata_thresholds={'GT_CONF':5})

def test_Genome_apply_vcf():

    indices=reference-sample_01
    assert indices[0]==4687
    assert reference.genome_coding_strand[reference.genome_index==indices[0]]=='t'
    assert sample_01.genome_coding_strand[reference.genome_index==indices[0]]=='c'
    assert sample_01.coverage[sample_01.genome_index==4687]==68
    assert sample_01.genome_sequence_metadata['GT_CONF'][sample_01.genome_index==4687][0]==pytest.approx(613.77)


    assert indices[1]==4725
    assert reference.genome_coding_strand[reference.genome_index==indices[1]]=='t'
    assert sample_01.genome_coding_strand[reference.genome_index==indices[1]]=='c'
    assert sample_01.coverage[sample_01.genome_index==4725]==68
    assert sample_01.genome_sequence_metadata['GT_CONF'][sample_01.genome_index==4725][0]==pytest.approx(613.77)

def test_Genome_list_variants_wrt():

    assert sample_01.list_variants_wrt(reference)==['4687t>c','4725t>c', '13333c>z','4730_indel','4735_indel']


def test_Genome_table_variants_wrt():

    # assert sample_01.table_variants_wrt(reference)==['4687t>c']
    foo=sample_01.table_variants_wrt(reference)
    pass


def test_Genome__complement():

    test_sequence=numpy.array(['a','c','t','g','z','x'])

    assert numpy.array_equal(reference._complement(test_sequence),numpy.array(['t','g','a','c','z','x']))

# use the subset argument to speed up the genome creation
h37rv=Genome(genbank_file="config/H37rV_v3.gbk",name="H37rV_v3",gene_subset=['katG','rpoB'])

def test_Genome_H37rV_katG():

    reversed_complement=''.join(h37rv.genome_sequence[h37rv.genome_feature_name=='katG'])[::-1]

    # below string taken from mycobrowser.epfl.ch and includes 100 bases upstream
    assert reversed_complement =='gtcatctactggggtctatgtcctgattgttcgatatccgacacttcgcgatcacatccgtgatcacagcccgataacaccaactcctggaaggaatgctgtgcccgagcaacacccacccattacagaaaccaccaccggagccgctagcaacggctgtcccgtcgtgggtcatatgaaataccccgtcgagggcggcggaaaccaggactggtggcccaaccggctcaatctgaaggtactgcaccaaaacccggccgtcgctgacccgatgggtgcggcgttcgactatgccgcggaggtcgcgaccatcgacgttgacgccctgacgcgggacatcgaggaagtgatgaccacctcgcagccgtggtggcccgccgactacggccactacgggccgctgtttatccggatggcgtggcacgctgccggcacctaccgcatccacgacggccgcggcggcgccgggggcggcatgcagcggttcgcgccgcttaacagctggcccgacaacgccagcttggacaaggcgcgccggctgctgtggccggtcaagaagaagtacggcaagaagctctcatgggcggacctgattgttttcgccggcaactgcgcgctggaatcgatgggcttcaagacgttcgggttcggcttcggccgggtcgaccagtgggagcccgatgaggtctattggggcaaggaagccacctggctcggcgatgagcgttacagcggtaagcgggatctggagaacccgctggccgcggtgcagatggggctgatctacgtgaacccggaggggccgaacggcaacccggaccccatggccgcggcggtcgacattcgcgagacgtttcggcgcatggccatgaacgacgtcgaaacagcggcgctgatcgtcggcggtcacactttcggtaagacccatggcgccggcccggccgatctggtcggccccgaacccgaggctgctccgctggagcagatgggcttgggctggaagagctcgtatggcaccggaaccggtaaggacgcgatcaccagcggcatcgaggtcgtatggacgaacaccccgacgaaatgggacaacagtttcctcgagatcctgtacggctacgagtgggagctgacgaagagccctgctggcgcttggcaatacaccgccaaggacggcgccggtgccggcaccatcccggacccgttcggcgggccagggcgctccccgacgatgctggccactgacctctcgctgcgggtggatccgatctatgagcggatcacgcgtcgctggctggaacaccccgaggaattggccgacgagttcgccaaggcctggtacaagctgatccaccgagacatgggtcccgttgcgagataccttgggccgctggtccccaagcagaccctgctgtggcaggatccggtccctgcggtcagccacgacctcgtcggcgaagccgagattgccagccttaagagccagatccgggcatcgggattgactgtctcacagctagtttcgaccgcatgggcggcggcgtcgtcgttccgtggtagcgacaagcgcggcggcgccaacggtggtcgcatccgcctgcagccacaagtcgggtgggaggtcaacgaccccgacggggatctgcgcaaggtcattcgcaccctggaagagatccaggagtcattcaactccgcggcgccggggaacatcaaagtgtccttcgccgacctcgtcgtgctcggtggctgtgccgccatagagaaagcagcaaaggcggctggccacaacatcacggtgcccttcaccccgggccgcacggatgcgtcgcaggaacaaaccgacgtggaatcctttgccgtgctggagcccaaggcagatggcttccgaaactacctcggaaagggcaacccgttgccggccgagtacatgctgctcgacaaggcgaacctgcttacgctcagtgcccctgagatgacggtgctggtaggtggcctgcgcgtcctcggcgcaaactacaagcgcttaccgctgggcgtgttcaccgaggcctccgagtcactgaccaacgacttcttcgtgaacctgctcgacatgggtatcacctgggagccctcgccagcagatgacgggacctaccagggcaaggatggcagtggcaaggtgaagtggaccggcagccgcgtggacctggtcttcgggtccaactcggagttgcgggcgcttgtcgaggtctatggcgccgatgacgcgcagccgaagttcgtgcaggacttcgtcgctgcctgggacaaggtgatgaacctcgacaggttcgacgtgcgctga'

    indices=h37rv.genome_index[h37rv.genome_feature_name=='katG']
    assert indices[0]==2153889
    assert indices[-1]==2156211

    positions=h37rv.genome_nucleotide_number[h37rv.genome_feature_name=='katG']

    # check there is no zero position
    assert numpy.sum(positions==0)==0
    assert positions[0]==2223
    assert positions[-1]==-100

def test_Genome_H37rV_rpoB():

    sequence=''.join(h37rv.genome_sequence[h37rv.genome_feature_name=='rpoB'])

    # below string taken from mycobrowser.epfl.ch and includes 100 bases upstream
    assert sequence =='cgccggccgaaaccgacaaaattatcgcggcgaacgggcccgtgggcaccgctcctctaagggctctcgttggtcgcatgaagtgctggaaggatgcatcttggcagattcccgccagagcaaaacagccgctagtcctagtccgagtcgcccgcaaagttcctcgaataactccgtacccggagcgccaaaccgggtctccttcgctaagctgcgcgaaccacttgaggttccgggactccttgacgtccagaccgattcgttcgagtggctgatcggttcgccgcgctggcgcgaatccgccgccgagcggggtgatgtcaacccagtgggtggcctggaagaggtgctctacgagctgtctccgatcgaggacttctccgggtcgatgtcgttgtcgttctctgaccctcgtttcgacgatgtcaaggcacccgtcgacgagtgcaaagacaaggacatgacgtacgcggctccactgttcgtcaccgccgagttcatcaacaacaacaccggtgagatcaagagtcagacggtgttcatgggtgacttcccgatgatgaccgagaagggcacgttcatcatcaacgggaccgagcgtgtggtggtcagccagctggtgcggtcgcccggggtgtacttcgacgagaccattgacaagtccaccgacaagacgctgcacagcgtcaaggtgatcccgagccgcggcgcgtggctcgagtttgacgtcgacaagcgcgacaccgtcggcgtgcgcatcgaccgcaaacgccggcaaccggtcaccgtgctgctcaaggcgctgggctggaccagcgagcagattgtcgagcggttcgggttctccgagatcatgcgatcgacgctggagaaggacaacaccgtcggcaccgacgaggcgctgttggacatctaccgcaagctgcgtccgggcgagcccccgaccaaagagtcagcgcagacgctgttggaaaacttgttcttcaaggagaagcgctacgacctggcccgcgtcggtcgctataaggtcaacaagaagctcgggctgcatgtcggcgagcccatcacgtcgtcgacgctgaccgaagaagacgtcgtggccaccatcgaatatctggtccgcttgcacgagggtcagaccacgatgaccgttccgggcggcgtcgaggtgccggtggaaaccgacgacatcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatggaccagaacaacccgctgtcggggttgacccacaagcgccgactgtcggcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgaggaggaccgccacgtggtggcacaggccaattcgccgatcgatgcggacggtcgcttcgtcgagccgcgcgtgctggtccgccgcaaggcgggcgaggtggagtacgtgccctcgtctgaggtggactacatggacgtctcgccccgccagatggtgtcggtggccaccgcgatgattcccttcctggagcacgacgacgccaaccgtgccctcatgggggcaaacatgcagcgccaggcggtgccgctggtccgtagcgaggccccgctggtgggcaccgggatggagctgcgcgcggcgatcgacgccggcgacgtcgtcgtcgccgaagaaagcggcgtcatcgaggaggtgtcggccgactacatcactgtgatgcacgacaacggcacccggcgtacctaccggatgcgcaagtttgcccggtccaaccacggcacttgcgccaaccagtgccccatcgtggacgcgggcgaccgagtcgaggccggtcaggtgatcgccgacggtccctgtactgacgacggcgagatggcgctgggcaagaacctgctggtggccatcatgccgtgggagggccacaactacgaggacgcgatcatcctgtccaaccgcctggtcgaagaggacgtgctcacctcgatccacatcgaggagcatgagatcgatgctcgcgacaccaagctgggtgcggaggagatcacccgcgacatcccgaacatctccgacgaggtgctcgccgacctggatgagcggggcatcgtgcgcatcggtgccgaggttcgcgacggggacatcctggtcggcaaggtcaccccgaagggtgagaccgagctgacgccggaggagcggctgctgcgtgccatcttcggtgagaaggcccgcgaggtgcgcgacacttcgctgaaggtgccgcacggcgaatccggcaaggtgatcggcattcgggtgttttcccgcgaggacgaggacgagttgccggccggtgtcaacgagctggtgcgtgtgtatgtggctcagaaacgcaagatctccgacggtgacaagctggccggccggcacggcaacaagggcgtgatcggcaagatcctgccggttgaggacatgccgttccttgccgacggcaccccggtggacattattttgaacacccacggcgtgccgcgacggatgaacatcggccagattttggagacccacctgggttggtgtgcccacagcggctggaaggtcgacgccgccaagggggttccggactgggccgccaggctgcccgacgaactgctcgaggcgcagccgaacgccattgtgtcgacgccggtgttcgacggcgcccaggaggccgagctgcagggcctgttgtcgtgcacgctgcccaaccgcgacggtgacgtgctggtcgacgccgacggcaaggccatgctcttcgacgggcgcagcggcgagccgttcccgtacccggtcacggttggctacatgtacatcatgaagctgcaccacctggtggacgacaagatccacgcccgctccaccgggccgtactcgatgatcacccagcagccgctgggcggtaaggcgcagttcggtggccagcggttcggggagatggagtgctgggccatgcaggcctacggtgctgcctacaccctgcaggagctgttgaccatcaagtccgatgacaccgtcggccgcgtcaaggtgtacgaggcgatcgtcaagggtgagaacatcccggagccgggcatccccgagtcgttcaaggtgctgctcaaagaactgcagtcgctgtgcctcaacgtcgaggtgctatcgagtgacggtgcggcgatcgaactgcgcgaaggtgaggacgaggacctggagcgggccgcggccaacctgggaatcaatctgtcccgcaacgaatccgcaagtgtcgaggatcttgcgtaa'

    indices=h37rv.genome_index[h37rv.genome_feature_name=='rpoB']
    assert indices[0]==759707
    assert indices[-1]==763325

    positions=h37rv.genome_nucleotide_number[h37rv.genome_feature_name=='rpoB']

    # check there is no zero position
    assert numpy.sum(positions==0)==0
    assert positions[0]==-100
    assert positions[-1]==3519


    # sample_04=copy.deepcopy(reference)
    # sample_04.apply_vcf_file(vcf_file=TEST_CASE_DIR+"04.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'])
    # (original_bases,indices,new_bases)=reference-sample_04
    # assert original_bases[0]=='c'
    # assert indices[0]==2155168
    # assert new_bases[0]=='g'
    # assert sample_04.coverage[sample_04.index==2155168]==53
    # assert sample_04.sequence_metadata['GT_CONF'][sample_04.index==2155168][0]==pytest.approx(500.23)
