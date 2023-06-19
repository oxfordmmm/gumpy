'''Tests relating to minor populations
'''
import pytest

import gumpy

#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def make_reproducable(l: [str]) -> [str]:
    '''Takes a list of tuples and sorted alphabetically so test case order doesn't matter

    Args:
        l ([str]): List to sort

    Returns:
        [str]: Sorted list
    '''
    #Convert each element of the list to a str, mapping to original
    str_l = {str(x): x for x in l}
    #Sort the str of each element
    sorted_str_l = sorted(list(str_l.keys()))
    #Map back to original and return
    return [str_l[key] for key in sorted_str_l]

@pytest.mark.slow
def test_minor_failures():
    '''Tests which are expected to fail/do nothing
    '''

    #VCF with minor populations, but no indices given
    vcf = gumpy.VCFFile(
                        'tests/test-cases/minor-populations.vcf',
                        ignore_filter=True, 
                        minor_population_indices={}
        )

    assert vcf.minor_populations == []

    ref = gumpy.Genome("config/NC_000962.3.gbk.gz", is_reference=True)

    sample = ref + vcf

    assert sample.minor_populations == []
    assert sample.minority_populations_GARC() == []

    #VCF which doesn't have minor populations
    #Only has 1 entry at 761155
    vcf = gumpy.VCFFile(
                        'tests/test-cases/03.vcf', 
                        minor_population_indices={761155}
    )
    assert vcf.minor_populations == []
    expected = ['VCF variant file, version 4.2', 'tests/test-cases/03.vcf', '1 records', 'FORMAT columns: COV, DP, GT, GT_CONF', "NC_000962.3 761155 c ('t',) 5.0 . GT:DP:COV:GT_CONF (1, 1):68:(0, 68):613.77"]
    assert [x.replace("\t", " ") for x in vcf.__repr__().split("\n") if x.replace("\t","")] == expected

    #Shouldn't pass the filter
    vcf = gumpy.VCFFile(
                        'tests/test-cases/03.vcf', 
                        minor_population_indices={761155},
                        format_fields_min_thresholds={'DP':100}
    )
    assert vcf.nucleotide_index.tolist() == []

    #Specifying minor indices outside of the genome shouldn't cause crashing
    vcf = gumpy.VCFFile(
                        'tests/test-cases/minor-populations.vcf',
                        ignore_filter=True,
                        minor_population_indices={9999999999999999999}
    )
    assert vcf.minor_populations == []
    sample = ref + vcf
    assert sample.minor_populations == []
    assert sample.minority_populations_GARC() == []

    #Expected exceptions
    with pytest.raises(Exception):
        vcf = gumpy.VCFFile(
                            'tests/test-cases/minor-populations.vcf',
                            ignore_filter=True,
                            minor_population_indices=[None]
        )
    with pytest.raises(Exception):
        vcf = gumpy.VCFFile(
                            'tests/test-cases/minor-populations.vcf',
                            ignore_filter=True,
                            minor_population_indices={None}
        )
    with pytest.raises(Exception):
        vcf = gumpy.VCFFile(
                            'tests/test-cases/minor-populations.vcf',
                            ignore_filter=True,
                            minor_population_indices='1,2,3'
        )
    with pytest.raises(Exception):
        #Giving it a file which isnt vcf
        vcf = gumpy.VCFFile(
                            'tests/test-cases/01.fasta',
                            ignore_filter=True,
                            minor_population_indices={1,2,3}
        )
    
    



@pytest.mark.slow
def test_get_minors():
    '''Tests of expected success for minor populations. Uses TB, FIXME if too slow
    '''
    vcf = gumpy.VCFFile(
                        'tests/test-cases/minor-populations.vcf',
                        ignore_filter=True, 
                        minor_population_indices={
                                                    7569, 7570, 7571, 7581, 7582, 7583, 7585,
                                                    7572, 7573, 7574
                        }
        )

    #Check for expected minority populations
    expected_minor_populations = [
                                    (7585, 'snp', ('g', 'g'), 15, 0.15), #Ref call
                                    (7582, 'del', 'ac', 10, 0.098), #Del
                                    (7581, 'snp', ('g', 't'), 10, 0.098), #SNP
                                    (7569, 'ins', 'gt', 10, 0.098), #Ins
                                    (7571, 'snp', ('g', 't'), 10, 0.098), #Synon SNP

                                    #Change all bases of a codon
                                    (7572, 'snp', ('t', 'a'), 2, 0.02), #SNP
                                    (7573, 'snp', ('c', 'c'), 10, 0.1), #Ref call
                                    (7574, 'snp', ('g', 't'), 10, 0.1), #SNP
        ]
    assert make_reproducable(vcf.minor_populations) == make_reproducable(expected_minor_populations)

    ref = gumpy.Genome("config/NC_000962.3.gbk.gz", is_reference=True, verbose=True)
    sample = ref + vcf

    #Checks for genome level minority populations
    assert sample.minority_populations_GARC() == sorted([
                                                        '7569_ins_gt:10', '7571g>t:10', '7581g>t:10', '7582_del_ac:10', '7585c>g:15',
                                                        '7572t>a:2', '7573t>c:10', '7574g>t:10'
                                                ])
    assert sample.minority_populations_GARC(reference=ref) == sorted([
                                                                    '7569_ins_gt:10', '7571g>t:10', '7581g>t:10', '7582_del_ac:10',
                                                                    '7585g>g:15', '7572t>a:2', '7573c>c:10', '7574g>t:10'
                                                                ])

    assert sample.minority_populations_GARC(interpretation='percentage') == sorted([
                                                                                    '7569_ins_gt:0.098', '7571g>t:0.098',
                                                                                    '7581g>t:0.098', '7582_del_ac:0.098', 
                                                                                    '7585c>g:0.15', '7572t>a:0.02', '7573t>c:0.1', 
                                                                                    '7574g>t:0.1'
                                                                            ])
    assert sample.minority_populations_GARC(interpretation='percentage', reference=ref) == sorted([
                                                                                                    '7569_ins_gt:0.098', '7571g>t:0.098',
                                                                                                    '7581g>t:0.098', '7582_del_ac:0.098', 
                                                                                                    '7585g>g:0.15', '7572t>a:0.02', 
                                                                                                    '7573c>c:0.1', '7574g>t:0.1'
                                                                                            ])

    #Check for gene level minority populations
    gyrA = sample.build_gene("gyrA")
    gyrA_ref = ref.build_gene("gyrA")

    assert gyrA.minority_populations_GARC() == sorted([
                                                        'D94Z:10', 'T95T:15',
                                                        'A90Z:10',
                                                        'L91Z:10'
                                                ])
    assert gyrA.minority_populations_GARC(reference=gyrA_ref) == sorted([
                                                                        'D94Z:10',
                                                                        'A90Z:10',
                                                                        'S91Z:10',
                                                                        'S95R:15'
                                                                    ])

    assert gyrA.minority_populations_GARC(interpretation='percentage') == sorted([
                                                                                'D94Z:0.098', 'T95T:0.15',
                                                                                'A90Z:0.098',
                                                                                'L91Z:0.1'
                                                                            ])
    assert gyrA.minority_populations_GARC(interpretation='percentage', reference=gyrA_ref) == sorted([
                                                                                                    'D94Z:0.098',
                                                                                                    'A90Z:0.098',
                                                                                                    'S91Z:0.1',
                                                                                                    'S95R:0.15'
                                                                                                ])
    
    #Checking for non-coding versions too (but hack around to make it so gyrA doesn't code)
    sample.genes['gyrA']['codes_protein'] = False
    gyrA = sample.build_gene("gyrA")
    assert sorted(gyrA.minority_populations_GARC()) == sorted([
                                                        '281_del_ac:10', '268_ins_gt:10',
                                                        'g270t:10', 'g280t:10', 'c284g:15',
                                                        't271a:2', 't272c:10', 'g273t:10'

                                                ])
    assert sorted(gyrA.minority_populations_GARC(reference=gyrA_ref)) == sorted([
                                                        '281_del_ac:10', '268_ins_gt:10',
                                                        'g270t:10', 'g280t:10',
                                                        't271a:2', 'g273t:10',
                                                        'c272c:10', 'g284g:15'
                                                ])    
    #Convert back to avoid problems
    sample.genes['gyrA']['codes_protein'] = True
    gyrA = sample.build_gene("gyrA")

    #Check GenomeDifference
    diff = ref - sample

    assert diff.minor_populations() == sample.minority_populations_GARC(reference=ref)
    assert diff.minor_populations(interpretation='percentage') == sample.minority_populations_GARC(reference=ref, interpretation='percentage')

    #Check GeneDifference
    geneDiff = gyrA_ref - gyrA
    assert geneDiff.minor_populations() == gyrA.minority_populations_GARC(reference=gyrA_ref)
    assert geneDiff.minor_populations(interpretation='percentage') == gyrA.minority_populations_GARC(reference=gyrA_ref, interpretation='percentage')

    #Should complain as the sample already has minor populations
    with pytest.raises(Exception):
        sample + vcf

    
    #Checking revcomp minors - use a different VCF with katG mutations for this
    #(tried hacking to make gyrA revcomp but this causes other issues as gene indices are already assigned)

    vcf = gumpy.VCFFile(
                        'tests/test-cases/minor-populations-revcomp.vcf',
                        ignore_filter=True, 
                        minor_population_indices=set(range(2154395, 2154410))
        )
    sample = ref + vcf

    katG = sample.build_gene("katG")
    assert sorted(katG.minority_populations_GARC()) == sorted([
                                                                'T572T:25', '1710_del_cc:25'
                                                            ])

    #Similarly, edge case of revcomp, non-coding (so hack katG to be non-coding)
    sample.genes['katG']['codes_protein'] = False
    katG = sample.build_gene("katG")
    assert sorted(katG.minority_populations_GARC()) == sorted([
                                                            'c1715a:25', '1710_del_cc:25'
                                                        ])
    