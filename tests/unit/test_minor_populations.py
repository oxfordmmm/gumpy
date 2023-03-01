import gumpy
import numpy
import pytest
import pickle

#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def test_get_minors():
    vcf = gumpy.VCFFile('tests/test-cases/minor-populations.vcf', ignore_filter=True, minor_population_indices={7569, 7570, 7571, 7581, 7582, 7583, 7585})

    ref = gumpy.Genome("config/NC_000962.3.gbk.gz", show_progress_bar=True, is_reference=True)
    sample = ref + vcf
    print()
    print("-----------------------")
    print()

    print(">",sample.minor_populations_GARC())
    print(">>",sample.minor_populations_GARC(interpretation='percentage'))


    # gyrA = sample.build_gene("gyrA")    

    assert False