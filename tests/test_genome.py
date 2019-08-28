import pytest, numpy, copy

from pathlib import Path

from gumpy import Genome

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

    assert len(reference.gene_names)==8

    # try a normal gene
    mask=reference.genome_feature_name=="M2"

    sequence=reference.genome_coding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tta"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="tag"

    sequence=reference.genome_noncoding_strand[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="aat"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="atc"

    sequence=reference.genome_sequence[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tta"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="tag"

reference2=Genome(fasta_file="config/NC_004148.2.fasta.gz",name="HMPV")

def test_Genome_instantiate_fasta():

    # check that the M. tuberculosis H37rV genome is the right length
    assert reference2.genome_length==13335

    # check the species is stored correctly
    assert reference2.organism=='Human metapneumovirus'

    # check that the sequence starts and ends as we expect
    assert reference2.genome_coding_strand[0]=='a'
    assert reference2.genome_coding_strand[-1]=='t'


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

# def test_Genome_contains_gene():
#
#     assert reference.contains_gene("rpoB")==True
#     assert reference.contains_gene("rpoD")==False
#     assert reference.contains_gene(5)==False
#     assert reference.contains_gene("")==False
#     assert reference.contains_gene("rpoBC")==False
#     assert reference.contains_gene("RPOB")==False
#
# def test_Genome_at_index():
#
#     assert reference.at_index(1471845)==('rrs', 'PROM')
#     assert reference.at_index(1471846)==('rrs', 'RNA')
#     assert reference.at_index(759807)==('rpoB', 'GENE')
#     assert reference.at_index(759806)==('rpoB', 'PROM')
#     assert reference.at_index(759707)==('rpoB', 'PROM')
#     assert reference.at_index(759706) is None
#     assert reference.at_index(763325)==('rpoB', 'GENE')
#     assert reference.at_index(763326)==('rpoC', 'PROM')
#
def test_Genome_calculate_snp_distance():

    sample=copy.deepcopy(reference)

    sample.genome_coding_strand[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based
    assert sample.snp_distance(reference)==1

    # reverse the change
    sample.genome_coding_strand[2]='g'
    assert sample.snp_distance(reference)==0

    sample.genome_coding_strand[2]='t'
    sample.genome_coding_strand[3]='t'
    assert sample.snp_distance(reference)==2

def test_Genome_apply_vcf():

    sample_01=copy.deepcopy(reference)
    sample_01.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF','GT_CONF_PERCENTILE'])
    indices=reference-sample_01
    assert indices[0]==4687
    assert reference.genome_coding_strand[reference.genome_index==indices[0]]=='t'
    assert sample_01.genome_coding_strand[reference.genome_index==indices[0]]=='c'
    assert sample_01.coverage[sample_01.genome_index==4687]==68
    assert sample_01.genome_sequence_metadata['GT_CONF'][sample_01.genome_index==4687][0]==pytest.approx(613.77)

    # sample_04=copy.deepcopy(reference)
    # sample_04.apply_vcf_file(vcf_file=TEST_CASE_DIR+"04.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'])
    # (original_bases,indices,new_bases)=reference-sample_04
    # assert original_bases[0]=='c'
    # assert indices[0]==2155168
    # assert new_bases[0]=='g'
    # assert sample_04.coverage[sample_04.index==2155168]==53
    # assert sample_04.sequence_metadata['GT_CONF'][sample_04.index==2155168][0]==pytest.approx(500.23)
