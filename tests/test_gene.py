import pytest, copy, numpy
from piezo import Gene

test_nucleotide_string="aaatttcccggg"
test_promoter_string="atcg"

# reference=Gene(   coding_nucleotides=test_nucleotide_string,\
#                         promoter_nucleotides=test_promoter_string,\
#                         gene_name="test",\
#                         first_nucleotide_index=100,\
#                         element_type='GENE',
#                         reverse=False )
#
# def test_Gene_instantiation():
#
#     assert reference.number_coding_nucleotides==len(test_nucleotide_string)
#     assert reference.number_promoter_nucleotides==len(test_promoter_string)
#     assert reference.gene_name=="test"
#     assert reference.reverse==False
#
# def test_Gene_hom_snp():
#
#     sample=copy.deepcopy(reference)
#     sample.apply_variant(position=100,ref_bases='a',alt_bases='c',coverage=20,model_score=30,model_percentile=40)
#     sample.update()
#     assert sample.coding_nucleotides_string==('caatttcccggg','caatttcccggg')
#     print(sample.coding_coverage)
#     assert numpy.array_equal(sample.coding_coverage,(numpy.array([20,0,0,0,0,0,0,0,0,0,0,0]),numpy.array([20,0,0,0,0,0,0,0,0,0,0,0])))
