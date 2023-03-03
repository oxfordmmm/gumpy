import gumpy
import numpy
import pytest
import pickle

#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

@pytest.mark.slow
def test_get_minors():
    vcf = gumpy.VCFFile('tests/test-cases/minor-populations.vcf', ignore_filter=True, minor_population_indices={7569, 7570, 7571, 7581, 7582, 7583, 7585})

    #TODO: Add ins + synonymous + 2 mutations in same codon + expected failures + expected no minor populations
    #Check for expected minority populations
    expected_minor_populations = [
                                    (7585, 'snp', ('g', 'g'), 15, 0.15), 
                                    (7582, 'del', 'ac', 10, 0.098), 
                                    (7581, 'snp', ('g', 't'), 10, 0.098)
        ]
    assert vcf.minor_populations == expected_minor_populations

    ref = gumpy.Genome("config/NC_000962.3.gbk.gz", is_reference=True)
    sample = ref + vcf

    #Checks for genome level minority populations
    assert sample.minority_populations_GARC() == ['7585c>g:15', '7582_del_ac:10', '7581g>t:10']
    assert sample.minority_populations_GARC(reference=ref) == ['7585g>g:15', '7582_del_ac:10', '7581g>t:10']

    assert sample.minority_populations_GARC(interpretation='percentage') == ['7585c>g:0.15', '7582_del_ac:0.098', '7581g>t:0.098']
    assert sample.minority_populations_GARC(interpretation='percentage', reference=ref) == ['7585g>g:0.15', '7582_del_ac:0.098', '7581g>t:0.098']

    #Check for gene level minority populations
    gyrA = sample.build_gene("gyrA")
    gyrA_ref = ref.build_gene("gyrA")
    assert gyrA.minority_populations_GARC() == ['gyrA@281_del_ac:10', 'gyrA@D94Y:10', 'gyrA@T95S:15']
    assert gyrA.minority_populations_GARC(reference=gyrA_ref) == ['gyrA@281_del_ac:10', 'gyrA@D94Y:10']

    assert gyrA.minority_populations_GARC(interpretation='percentage') == ['gyrA@281_del_ac:0.098', 'gyrA@D94Y:0.098', 'gyrA@T95S:0.15']
    assert gyrA.minority_populations_GARC(interpretation='percentage', reference=gyrA_ref) == ['gyrA@281_del_ac:0.098', 'gyrA@D94Y:0.098']

    #Check GenomeDifference
    diff = sample - ref

    assert diff.minor_populations() == sample.minority_populations_GARC(reference=ref)
    assert diff.minor_populations(interpretation='percentage') == sample.minority_populations_GARC(reference=ref, interpretation='percentage')

    #Check GeneDifference
    geneDiff = gyrA - gyrA_ref
    assert geneDiff.minor_populations() == gyrA.minority_populations_GARC(reference=gyrA_ref)
    assert geneDiff.minor_populations(interpretation='percentage') == gyrA.minority_populations_GARC(reference=gyrA_ref, interpretation='percentage')