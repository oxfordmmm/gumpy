import pytest, copy, numpy
from gumpy.genome import Genome

test_nucleotide_string="aaatttcccggg"
test_promoter_string="atcg"

reference=Genome(genbank_file="config/NC_004148.2.gbk",name="HMPV")

test_gene=reference.genes["F"]

def test_Gene_instantiation():

    assert test_gene.total_number_nucleotides==1720
    assert test_gene.gene_name=="F"
    assert test_gene.on_noncoding_strand==False
    assert test_gene.numbering[0]==-100
    assert test_gene.numbering[-1]==540
    assert test_gene.numbering[-2]==540
    assert test_gene.index[test_gene.position==1]==3067
    assert test_gene.index[test_gene.position==2]==3068
    assert test_gene.index[test_gene.position==-1]==3066
    assert test_gene.sequence[0]=='a'
    assert test_gene.sequence[-1]=='g'
    assert test_gene.sequence[-2]=='a'

    F_gene="MSWKVVIIFSLLITPQHGLKESYLEESCSTITEGYLSVLRTGWYTNVFTLEVGDVENLTCSDGPSLIKTELDLTKSALRELKTVSADQLAREEQIENPRQSRFVLGAIALGVATAAAVTAGVAIAKTIRLESEVTAIKNALKTTNEAVSTLGNGVRVLATAVRELKDFVSKNLTRAINKNKCDIDDLKMAVSFSQFNRRFLNVVRQFSDNAGITPAISLDLMTDAELARAVSNMPTSAGQIKLMLENRAMVRRKGFGILIGVYGSSVIYMVQLPIFGVIDTPCWIVKAAPSCSGKKGNYACLLREDQGWYCQNAGSTVYYPNEKDCETRGDHVFCDTAAGINVAEQSKECNINISTTNYPCKVSTGRHPISMVALSPLGALVACYKGVSCSIGSNRVGIIKQLNKGCSYITNQDADTVTIDNTVYQLSKVEGEQHVIKGRPVSSSFDPIKFPEDQFNVALDQVFENIENSQALVDQSNRILSSAEKGNTGFIIVIILIAVLGSSMILVSIFIIIKKTKKPTGAPPELSGVTNNGFIPHS!"

    assert F_gene==''.join(test_gene.amino_acid_sequence)


# def test_Gene___repr__():
#
#     expected_output="F gene\n1720 nucleotides, codes for protein\nactac...taaaa\n-100 -99 -98 -97 -96 ...-5 -4 -3 -2 -1\nMSWKV...IPHS!\n1 2 3 4 5 ...536 537 538 539 540 \n"
#
#     assert test_gene.__repr__()==expected_output

new_gene=copy.deepcopy(test_gene)
new_gene.sequence[new_gene.position==1]='t'
new_gene._translate_sequence()

def test_Gene_list_mutations_wrt_1():

    assert new_gene.list_mutations_wrt(test_gene)==['M1L']

def test_Gene___sub__1():

    assert new_gene-test_gene==[1]

new_gene2=copy.deepcopy(test_gene)
new_gene2.sequence[new_gene2.position==1]='c'
new_gene2.sequence[new_gene2.position==4]='c'
new_gene2._translate_sequence()

def test_Gene_list_mutations_wrt_2():

    assert new_gene2.list_mutations_wrt(test_gene)==['M1L','S2P']

def test_Gene___sub__2():

    assert list(new_gene2-test_gene)==[1,2]


# mutate it away, and then back again
new_gene3=copy.deepcopy(test_gene)
new_gene3.sequence[new_gene3.position==1]='t'
new_gene3._translate_sequence()
new_gene3.sequence[new_gene3.position==1]='a'
new_gene3._translate_sequence()

def test_Gene_list_mutations_wrt_3():

    assert new_gene3.list_mutations_wrt(test_gene) is None

def test_Gene_valid_variant():

    # badly formed variant
    with pytest.raises(Exception):
        assert test_gene.valid_variant("1a>g")

    assert test_gene.valid_variant("a1g")
    assert test_gene.valid_variant("t2g")
    assert test_gene.valid_variant("g3a")
    assert test_gene.valid_variant("g3z")
    assert test_gene.valid_variant("g3?")
    assert test_gene.valid_variant("a-100g")
    assert test_gene.valid_variant("g1620a")

    # out of range position
    with pytest.raises(Exception):
        assert test_gene.valid_variant("a1621g")
    with pytest.raises(Exception):
        assert test_gene.valid_variant("a-101g")

    # badly formed reference base
    with pytest.raises(Exception):
        assert test_gene.valid_variant("A1g")
    with pytest.raises(Exception):
        assert test_gene.valid_variant("p2g")
    with pytest.raises(Exception):
        assert test_gene.valid_variant("z2g")

    # badly formed positions
    with pytest.raises(Exception):
        assert test_gene.valid_variant("zag")
    with pytest.raises(Exception):
        assert test_gene.valid_variant("z1o1g")
