import numpy, gumpy, pytest, math


#As BioPython thinks that the locus line of the TEST-RNA.gbk is malformed, it gives a warning
#So ignore it to stop failing tests...
pytestmark = pytest.mark.filterwarnings("ignore")

def test_instanciate_genome_rna():
    genome = gumpy.Genome("config/TEST-RNA.gbk")
    g2 = gumpy.Genome("config/TEST-RNA.gbk", multithreaded=True)
    #Check that multithreading gives the same result.
    #If not testing on Linux the multithreaded code will never be touched...
    assert genome == g2

    #Testing generic attributes such as name and length
    assert len(genome) == 99
    assert genome.name == "TEST_RNA"
    assert genome.id == "TEST_RNA.1"
    assert list(genome.genes_lookup.keys()) == list("ABC")

    #Testing annotations parsed correctly
    assert genome.annotations['organism'] == "TEST_RNA_ORGANISM"
    assert genome.annotations['source'] == "TEST_RNA_SOURCE"
    assert genome.annotations['references'] == [{
                                                "location": [{
                                                     "_start": 0,
                                                     "_end": 99
                                                    }],
                                                "authors": "Test,1., Test,2., Test,3.",
                                                "consrtm": "",
                                                "title": "Test title for a reference for an RNA strand",
                                                "journal": "Test journal for an RNA strand",
                                                "medline_id": "",
                                                "pubmed_id": "1",
                                                "comment": ""
                                                 }]
    #Testing sequence was parsed correctly
    assert numpy.all(
        genome.nucleotide_sequence == numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        ))
    assert numpy.all(
        genome.nucleotide_index == numpy.array(range(1,100)))

    #Testing all arrays were setup correctly
    #Stacked arrays are changed during promoter assignment, so should test promoter assignment too
    original_stacked_gene_name = numpy.array([
        ['', '', '', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', '', '', ''],
        ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    ])
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    assert numpy.all(
        genome.stacked_gene_name == full_gene_name)
    #No rev-comp as RNA only has a single strand
    assert numpy.all(
        genome.stacked_is_reverse_complement == numpy.array([[False for i in range(99)], [False for i in range(99)]])
    )
    assert numpy.all(
        genome.stacked_nucleotide_index == numpy.array([[i for i in range(1, 100)], [i for i in range(1, 100)]])
    )
    assert numpy.all(
        genome.stacked_is_promoter == numpy.array([
            [True, True, True]+[False for i in range(3, 60)]+[True for i in range(60, 90)]+[False, False, False, False, False, False, True, True, True],
            [False for i in range(99)]
        ])
    )
    assert numpy.all(
        genome.stacked_nucleotide_number == numpy.array([
            [-3, -2, -1] + [i for i in range(1, 28)] + [0 for i in range(28, 58)] + [i for i in range(-30, 0)] + [i for i in range(1, 7)] + [-6, -5, -4],
            [0 for i in range(27)] + [i for i in range(1, 34)] + [0 for i in range(39)]
        ])
    )
    assert numpy.all(
        genome.stacked_is_cds == (original_stacked_gene_name != "")
    )

    #Testing all genes instanciated (not correct genes)
    assert list(genome.genes.keys()) == ["A", "B", "C"]

    assert genome.genes["A"] == gumpy.Gene(
        name="A",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "A"],
        index=genome.nucleotide_index[full_gene_name[0] == "A"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "A"],
        is_cds=genome.stacked_is_cds[full_gene_name == "A"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "A"],
        is_indel=genome.is_indel[full_gene_name[0] == "A"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )
    assert genome.genes["B"] == gumpy.Gene(
        name="B",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[1] == "B"],
        index=genome.nucleotide_index[full_gene_name[1] == "B"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "B"],
        is_cds=genome.stacked_is_cds[full_gene_name == "B"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "B"],
        is_indel=genome.is_indel[full_gene_name[1] == "B"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )

    assert genome.genes["C"] == gumpy.Gene(
        name="C",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "C"],
        index=genome.nucleotide_index[full_gene_name[0] == "C"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "C"],
        is_cds=genome.stacked_is_cds[full_gene_name == "C"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "C"],
        is_indel=genome.is_indel[full_gene_name[0] == "C"],
        indel_length=[0 for i in range(36)],
        feature_type="GENE"
    )

def test_instanciate_genes_rna():
    genome = gumpy.Genome("config/TEST-RNA.gbk")

    #If the previous checks pass, this should be equal to a freshly instanciated Gene
    gene = genome.genes["A"]

    #Ground truth values for all genes - used to extract values as required
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'A', 'A', 'A'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    #Nucleotides for the genome
    nucleotide_sequence = numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        )
    #Indices for nucleotides
    nucleotide_index = numpy.array(range(1,100))

    #Test generic attributes which are passed into the constructor
    assert gene.name == "A"
    assert gene.feature_type == "GENE"
    assert numpy.all(gene.nucleotide_sequence == nucleotide_sequence[full_gene_name[0] == "A"])
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))+[-6, -5, -4]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]+[False, False, False]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]+[True, True, True]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(33)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(33)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

def test_instanciate_genome_dna():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    #Testing generic attributes such as name and length
    assert len(genome) == 99
    assert genome.name == "TEST_DNA"
    assert genome.id == "TEST_DNA.1"
    assert list(genome.genes_lookup.keys()) == list("ABC")

    #Testing annotations parsed correctly
    assert genome.annotations['organism'] == "TEST_DNA_ORGANISM"
    assert genome.annotations['source'] == "TEST_DNA_SOURCE"
    assert genome.annotations['references'] == [{
                                                "location": [{
                                                     "_start": 0,
                                                     "_end": 99
                                                    }],
                                                "authors": "Test,1., Test,2., Test,3.",
                                                "consrtm": "",
                                                "title": "Test title for a reference for an DNA strand",
                                                "journal": "Test journal for an DNA strand",
                                                "medline_id": "",
                                                "pubmed_id": "1",
                                                "comment": ""
                                                 }]
    #Testing sequence was parsed correctly
    assert numpy.all(
        genome.nucleotide_sequence == numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        ))
    assert numpy.all(
        genome.nucleotide_index == numpy.array(range(1,100)))

    #Testing all arrays were setup correctly
    #Stacked arrays are changed during promoter assignment, so should test promoter assignment too
    original_stacked_gene_name = numpy.array([
        ['', '', '', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', '', '', ''],
        ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
    ])
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    assert numpy.all(
        genome.stacked_gene_name == full_gene_name)
    #No rev-comp as DNA only has a single strand
    assert numpy.all(
        genome.stacked_is_reverse_complement == numpy.array([[False for i in range(90)]+[True for i in range(9)], [False for i in range(99)]])
    )
    assert numpy.all(
        genome.stacked_nucleotide_index == numpy.array([[i for i in range(1, 100)], [i for i in range(1, 100)]])
    )
    assert numpy.all(
        genome.stacked_is_promoter == numpy.array([
            [True, True, True]+[False for i in range(3, 90)]+[False, False, False, False, False, False, True, True, True],
            [False for i in range(99)]
        ])
    )
    assert numpy.all(
        genome.stacked_nucleotide_number == numpy.array([
            [-3, -2, -1] + [i for i in range(1, 28)] + [0 for i in range(28, 88)] + [i for i in range(1, 7)][::-1] + [-1, -2, -3],
            [0 for i in range(27)] + [i for i in range(1, 34)] + [0 for i in range(39)]
        ])
    )
    assert numpy.all(
        genome.stacked_is_cds == (original_stacked_gene_name != "")
    )

    #Testing all genes instanciated (not correct genes)
    assert list(genome.genes.keys()) == ["A", "B", "C"]

    assert genome.genes["A"] == gumpy.Gene(
        name="A",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "A"],
        index=genome.nucleotide_index[full_gene_name[0] == "A"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "A"],
        is_cds=genome.stacked_is_cds[full_gene_name == "A"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "A"],
        is_indel=genome.is_indel[full_gene_name[0] == "A"],
        indel_length=[0 for i in range(30)],
        feature_type="GENE"
    )
    assert genome.genes["B"] == gumpy.Gene(
        name="B",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[1] == "B"],
        index=genome.nucleotide_index[full_gene_name[1] == "B"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "B"],
        is_cds=genome.stacked_is_cds[full_gene_name == "B"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "B"],
        is_indel=genome.is_indel[full_gene_name[1] == "B"],
        indel_length=[0 for i in range(33)],
        feature_type="GENE"
    )
    assert genome.genes["C"] == gumpy.Gene(
        name="C",
        nucleotide_sequence=genome.nucleotide_sequence[full_gene_name[0] == "C"],
        index=genome.nucleotide_index[full_gene_name[0] == "C"],
        nucleotide_number=genome.stacked_nucleotide_number[full_gene_name == "C"],
        is_cds=genome.stacked_is_cds[full_gene_name == "C"],
        is_promoter=genome.stacked_is_promoter[full_gene_name == "C"],
        is_indel=genome.is_indel[full_gene_name[0] == "C"],
        indel_length=[0 for i in range(9)],
        feature_type="GENE",
        reverse_complement=True
    )

def test_instanciate_genes_dna():
    genome = gumpy.Genome("config/TEST-DNA.gbk")

    #If the previous checks pass, this should be equal to a freshly instanciated Gene
    gene = genome.genes["A"]

    #Ground truth values for all genes - used to extract values as required
    #Grown out genes with promoters
    full_gene_name = numpy.array([
                ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ])
    #Nucleotides for the genome
    nucleotide_sequence = numpy.array(
            list("aaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccc")
        )
    #Indices for nucleotides
    nucleotide_index = numpy.array(range(1,100))

    #Test generic attributes which are passed into the constructor
    assert gene.name == "A"
    assert gene.feature_type == "GENE"
    assert numpy.all(gene.nucleotide_sequence == nucleotide_sequence[full_gene_name[0] == "A"])
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "A"])
    assert numpy.all(gene.nucleotide_number == numpy.array([-3, -2, -1]+list(range(1, 28))))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(27)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(27)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(30)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(30)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 28)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,10)))
    assert numpy.all(gene.codons == numpy.array(["aaa", "aaa", "acc", "ccc", "ccc", "ccg", "ggg", "ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("KKTPPPGGG")))

    #Checks for a reverse complement gene
    complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z', 'o':'o', 'n':'n', 'r':'y', 'y':'r', 's':'w', 'w':'s'}
    gene = genome.genes["C"]

    #Test generic attributes which are passed into the constructor
    assert gene.name == "C"
    assert gene.feature_type == "GENE"
    assert gene.reverse_complement == True
    assert numpy.all(gene.nucleotide_sequence == [complementary_bases[n] for n in nucleotide_sequence[full_gene_name[0] == "C"]])
    assert numpy.all(gene.index == nucleotide_index[full_gene_name[0] == "C"][::-1])
    assert numpy.all(gene.nucleotide_number == numpy.array((list(range(1, 7))[::-1]+[-1, -2, -3])[::-1]))
    assert numpy.all(gene.is_cds == numpy.array([False, False, False]+[True for i in range(6)]))
    assert numpy.all(gene.is_promoter == numpy.array([True, True, True]+[False for i in range(6)]))
    assert numpy.all(gene.is_indel == numpy.array([False for i in range(9)]))
    assert numpy.all(gene.indel_length == numpy.array([0 for i in range(9)]))

    #Test attributes set during object instanciation
    assert numpy.all(gene.triplet_number == numpy.array([math.ceil(i/3) for i in range(1, 7)]))
    assert numpy.all(gene.amino_acid_number == numpy.array(range(1,3)))
    assert numpy.all(gene.codons == numpy.array(["ggg", "ggg"]))
    assert numpy.all(gene.amino_acid_sequence == numpy.array(list("GG")))

def test_instanciate_vcf():
    vcf = gumpy.VariantFile("tests/test-cases/TEST-DNA.vcf")
    #Testing some populated items for the object
    assert vcf.VCF_VERSION == (4, 2)
    assert vcf.contig_lengths == {"TEST_DNA": 99}
    assert vcf.formats == {
        "COV": {
            "type":"Integer",
            "description":"Number of reads on ref and alt alleles",
            "id": 1
        },
        "GT": {
            "type": "String",
            "description": "Genotype",
            "id": 2
        },
        "DP": {
            "type": "Integer",
            "description": "total kmer depth from gramtools",
            "id": 3
        },
        "GT_CONF": {
            "type": "Float",
            "description": "Genotype confidence. Difference in log likelihood of most likely and next most likely genotype",
            "id": 4
        }
    }
    assert len(vcf.records) == 6

    #Due to the dict structure here, several asserts are required
    changes = {
        1: ('g', [(0, '*'), (68, 'g')]),
        15: ('t', [(0, '*'), (68, 't')]),
        27: ('z', [(1, '*'), (99, 't'), (100, 'c')]),
        71: (numpy.array(['g', 'c', 'c']), [(0, '*'), (68, numpy.array(['g', 'c', 'c']))]),
        77: ('z', [(0, '*'), (48, numpy.array(['g', 't', 't'])), (20, 'g')]),
        89: ('x', [(0, '*'), (68, 'x')])
    }
    assert vcf.changes.keys() == changes.keys()
    for key in changes.keys():
        assert numpy.all(vcf.changes[key][0] == changes[key][0])
        for ((n_reads1, call1), (n_reads2, call2)) in zip(vcf.changes[key][1], changes[key][1]):
            assert n_reads1 == n_reads2
            assert numpy.all(call1 == call2)

    #Testing record objects

    #Features common to all record objects:
    for record in vcf.records:
        assert record.chrom == "TEST_DNA"
    #Pos
    assert vcf.records[0].pos == 2
    assert vcf.records[1].pos == 16
    assert vcf.records[2].pos == 28
    assert vcf.records[3].pos == 72
    assert vcf.records[4].pos == 78
    assert vcf.records[5].pos == 90

    #Ref
    assert vcf.records[0].ref == "A"
    assert vcf.records[1].ref == "C"
    assert vcf.records[2].ref == "G"
    assert vcf.records[3].ref == "T"
    assert vcf.records[4].ref == "T"
    assert vcf.records[5].ref == "A"

    #Alt
    assert numpy.all(vcf.records[0].alts == ("G", ))
    assert numpy.all(vcf.records[1].alts == ("T", ))
    assert numpy.all(vcf.records[2].alts == ("T", "C"))
    assert numpy.all(vcf.records[3].alts == ("GCC", ))
    assert numpy.all(vcf.records[4].alts == ("GTT", "G"))
    assert numpy.all(vcf.records[5].alts == None)

    #Qual
    assert vcf.records[0].qual == None
    assert vcf.records[1].qual == None
    assert vcf.records[2].qual == None
    assert vcf.records[3].qual == None
    assert vcf.records[4].qual == None
    assert vcf.records[5].qual == None

    #Filter
    assert vcf.records[0].filter == "PASS"
    assert vcf.records[1].filter == "."
    assert vcf.records[2].filter == "PASS"
    assert vcf.records[3].filter == "PASS"
    assert vcf.records[4].filter == "."
    assert vcf.records[5].filter == "."

    #Other fields
    assert vcf.records[0].info == {'KMER': 15}
    assert vcf.records[1].info == {'KMER': 15}
    assert vcf.records[2].info == {'KMER': 15}
    assert vcf.records[3].info == {'KMER': 15}
    assert vcf.records[4].info == {'KMER': 15}
    assert vcf.records[5].info == {'KMER': 15}

    print(vcf.records[0].values.keys())
    #GT
    gt = [record.values["GT"] for record in vcf.records]
    #None is given as a GT value for null values for alts
    assert numpy.all(gt == [(1, 1), (1, 1), (1, 2), (1, 1), (1, 1), (None, None)])
    #DP
    dp = [record.values["DP"] for record in vcf.records]
    assert numpy.all(dp == [68, 68, 200, 68, 68, 68])
    #COV
    cov = [record.values["COV"] for record in vcf.records]
    assert numpy.all(cov == [(0, 68), (0, 68), (1, 99, 100), (0, 68), (0, 48, 20), (0, 68)])
    #GT_CONF
    gt_conf = [record.values["GT_CONF"] for record in vcf.records]
    assert numpy.all(gt_conf == [613.77, 613.77, 613.77, 63.77, 63.77, 63.77])

    #Quick test for VCFRecord.__repr__()
    assert vcf.records[0].__repr__() == "TEST_DNA\t2\tA\t('G',)\tNone\tPASS\tGT:DP:COV:GT_CONF\t(1, 1):68:(0, 68):613.77\n"


# import pytest, numpy, copy

# from pathlib import Path

# from gumpy.genome import Genome 

# TEST_CASE_DIR = "tests/test-cases/"

# reference=Genome("config/NC_045512.2.gbk")

# def test_Genome_instantiate_genbank():

#     # check that the SARS-CoV-2 genome is the right length
#     assert reference.length==29903

#     assert reference.name=='NC_045512'

#     assert reference.id=='NC_045512.2'

#     # check the species is stored correctly
#     assert reference.annotations['organism']=='Severe acute respiratory syndrome coronavirus 2'

#     # check that the sequence starts and ends as we expect
#     assert reference.nucleotide_sequence[0]=='a'
#     assert reference.nucleotide_sequence[-1]=='a'

#     # check the number of genes
#     assert len(reference.genes_lookup.keys())==12

#     # try a normal gene
#     mask=(reference.stacked_gene_name=="S") & (reference.stacked_is_cds)
#     sequence=reference.nucleotide_sequence[numpy.any(mask,axis=0)]
#     first_codon="".join(i for i in sequence[:3])
#     assert first_codon=="atg"
#     last_codon="".join(i for i in sequence[-3:])
#     assert last_codon=="taa"

#     spike_sequence=' '.join(i for i in sequence)

#     assert spike_sequence=='atgtttgtttttcttgttttattgccactagtctctagtcagtgtgttaatcttacaaccagaactcaattaccccctgcatacactaattctttcacacgtggtgtttattaccctgacaaagttttcagatcctcagttttacattcaactcaggacttgttcttacctttcttttccaatgttacttggttccatgctatacatgtctctgggaccaatggtactaagaggtttgataaccctgtcctaccatttaatgatggtgtttattttgcttccactgagaagtctaacataataagaggctggatttttggtactactttagattcgaagacccagtccctacttattgttaataacgctactaatgttgttattaaagtctgtgaatttcaattttgtaatgatccatttttgggtgtttattaccacaaaaacaacaaaagttggatggaaagtgagttcagagtttattctagtgcgaataattgcacttttgaatatgtctctcagccttttcttatggaccttgaaggaaaacagggtaatttcaaaaatcttagggaatttgtgtttaagaatattgatggttattttaaaatatattctaagcacacgcctattaatttagtgcgtgatctccctcagggtttttcggctttagaaccattggtagatttgccaataggtattaacatcactaggtttcaaactttacttgctttacatagaagttatttgactcctggtgattcttcttcaggttggacagctggtgctgcagcttattatgtgggttatcttcaacctaggacttttctattaaaatataatgaaaatggaaccattacagatgctgtagactgtgcacttgaccctctctcagaaacaaagtgtacgttgaaatccttcactgtagaaaaaggaatctatcaaacttctaactttagagtccaaccaacagaatctattgttagatttcctaatattacaaacttgtgcccttttggtgaagtttttaacgccaccagatttgcatctgtttatgcttggaacaggaagagaatcagcaactgtgttgctgattattctgtcctatataattccgcatcattttccacttttaagtgttatggagtgtctcctactaaattaaatgatctctgctttactaatgtctatgcagattcatttgtaattagaggtgatgaagtcagacaaatcgctccagggcaaactggaaagattgctgattataattataaattaccagatgattttacaggctgcgttatagcttggaattctaacaatcttgattctaaggttggtggtaattataattacctgtatagattgtttaggaagtctaatctcaaaccttttgagagagatatttcaactgaaatctatcaggccggtagcacaccttgtaatggtgttgaaggttttaattgttactttcctttacaatcatatggtttccaacccactaatggtgttggttaccaaccatacagagtagtagtactttcttttgaacttctacatgcaccagcaactgtttgtggacctaaaaagtctactaatttggttaaaaacaaatgtgtcaatttcaacttcaatggtttaacaggcacaggtgttcttactgagtctaacaaaaagtttctgcctttccaacaatttggcagagacattgctgacactactgatgctgtccgtgatccacagacacttgagattcttgacattacaccatgttcttttggtggtgtcagtgttataacaccaggaacaaatacttctaaccaggttgctgttctttatcaggatgttaactgcacagaagtccctgttgctattcatgcagatcaacttactcctacttggcgtgtttattctacaggttctaatgtttttcaaacacgtgcaggctgtttaataggggctgaacatgtcaacaactcatatgagtgtgacatacccattggtgcaggtatatgcgctagttatcagactcagactaattctcctcggcgggcacgtagtgtagctagtcaatccatcattgcctacactatgtcacttggtgcagaaaattcagttgcttactctaataactctattgccatacccacaaattttactattagtgttaccacagaaattctaccagtgtctatgaccaagacatcagtagattgtacaatgtacatttgtggtgattcaactgaatgcagcaatcttttgttgcaatatggcagtttttgtacacaattaaaccgtgctttaactggaatagctgttgaacaagacaaaaacacccaagaagtttttgcacaagtcaaacaaatttacaaaacaccaccaattaaagattttggtggttttaatttttcacaaatattaccagatccatcaaaaccaagcaagaggtcatttattgaagatctacttttcaacaaagtgacacttgcagatgctggcttcatcaaacaatatggtgattgccttggtgatattgctgctagagacctcatttgtgcacaaaagtttaacggccttactgttttgccacctttgctcacagatgaaatgattgctcaatacacttctgcactgttagcgggtacaatcacttctggttggacctttggtgcaggtgctgcattacaaataccatttgctatgcaaatggcttataggtttaatggtattggagttacacagaatgttctctatgagaaccaaaaattgattgccaaccaatttaatagtgctattggcaaaattcaagactcactttcttccacagcaagtgcacttggaaaacttcaagatgtggtcaaccaaaatgcacaagctttaaacacgcttgttaaacaacttagctccaattttggtgcaatttcaagtgttttaaatgatatcctttcacgtcttgacaaagttgaggctgaagtgcaaattgataggttgatcacaggcagacttcaaagtttgcagacatatgtgactcaacaattaattagagctgcagaaatcagagcttctgctaatcttgctgctactaaaatgtcagagtgtgtacttggacaatcaaaaagagttgatttttgtggaaagggctatcatcttatgtccttccctcagtcagcacctcatggtgtagtcttcttgcatgtgacttatgtccctgcacaagaaaagaacttcacaactgctcctgccatttgtcatgatggaaaagcacactttcctcgtgaaggtgtctttgtttcaaatggcacacactggtttgtaacacaaaggaatttttatgaaccacaaatcattactacagacaacacatttgtgtctggtaactgtgatgttgtaataggaattgtcaacaacacagtttatgatcctttgcaacctgaattagactcattcaaggaggagttagataaatattttaagaatcatacatcaccagatgttgatttaggtgacatctctggcattaatgcttcagttgtaaacattcaaaaagaaattgaccgcctcaatgaggttgccaagaatttaaatgaatctctcatcgatctccaagaacttggaaagtatgagcagtatataaaatggccatggtacatttggctaggttttatagctggcttgattgccatagtaatggtgacaattatgctttgctgtatgaccagttgctgtagttgtctcaagggctgttgttcttgtggatcctgctgcaaatttgatgaagacgactctgagccagtgctcaaaggagtcaaattacattacacataa'



# # reference2=Genome(fasta_file="config/NC_004148.2.fasta.gz",name="HMPV")
# #
# # def test_Genome_instantiate_fasta():
# #
# #     # check that the M. tuberculosis H37rV genome is the right length
# #     assert reference2.genome_length==13335
# #
# #     # check the species is stored correctly
# #     assert reference2.organism=='Human metapneumovirus'
# #
# #     # check that the sequence starts and ends as we expect
# #     assert reference2.genome_coding_strand[0]=='a'
# #     assert reference2.genome_coding_strand[-1]=='t'
# #
# # def test_Genome_valid_gene_mutation_snps():
# #
# #     # correct protein SNPs
# #     assert reference.valid_gene_mutation("F@M1N")
# #     assert reference.valid_gene_mutation("F@M1Z")
# #     assert reference.valid_gene_mutation("F@M1O")
# #     assert reference.valid_gene_mutation("F@M1X")
# #     assert reference.valid_gene_mutation("F@S2?")
# #     assert reference.valid_gene_mutation("F@S2=")
# #     assert reference.valid_gene_mutation("F@S2!")
# #     assert reference.valid_gene_mutation("F@a-1c")
# #     assert reference.valid_gene_mutation("F@a-1o")
# #     assert reference.valid_gene_mutation("F@a-1z")
# #     assert reference.valid_gene_mutation("F@a-1x")
# #     assert reference.valid_gene_mutation("M2@T76L")
# #     assert reference.valid_gene_mutation("M2@*?")
# #     assert reference.valid_gene_mutation("M2@*Z")
# #     assert reference.valid_gene_mutation("M2@*X")
# #     assert reference.valid_gene_mutation("M2@*O")
# #     assert reference.valid_gene_mutation("M2@-*?")
# #     assert reference.valid_gene_mutation("M2@-*o")
# #     assert reference.valid_gene_mutation("M2@-*z")
# #     assert reference.valid_gene_mutation("M2@-*x")
# #
# #     # just badly formed
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("____")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("o_o_o_o_9")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("flkgjslkjg")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("o_i9k")
# #
# #     # genes not present
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("lkdfjlksdjf_P1N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("rpoB@P76L")
# #
# #     # incorrect reference amino acids
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@P1N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("M2@P76L")
# #
# #     # bad reference amino acids
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@;1N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@B1N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@81N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("M2@J76L")
# #
# #     # bad positions
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@PKN")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@P-2N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@P1000N")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@P:N")
# #
# #     # bad target amino acids
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@M1OO")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@PKB")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@PK;")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@PKJ")
# #
# #
# # def test_Genome_valid_gene_mutation_indels():
# #
# #     # correct INDELs with good grammar
# #     assert reference.valid_gene_mutation("F@1_indel")
# #     assert reference.valid_gene_mutation("F@1_ins")
# #     assert reference.valid_gene_mutation("F@1_del")
# #     assert reference.valid_gene_mutation("F@1_ins_3")
# #     assert reference.valid_gene_mutation("F@1_ins_ctga")
# #     assert reference.valid_gene_mutation("F@1_fs")
# #     assert reference.valid_gene_mutation("F@1_del_3")
# #     assert reference.valid_gene_mutation("F@-1_indel")
# #     assert reference.valid_gene_mutation("F@1_del_acgt")
# #
# #
# #     # bad grammar
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@1_indl")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@1_frameshift")
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@1_ins_ggaf")
# #
# #     # wrong ordering
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@indel_1")
# #
# #     # incorrect gene
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F1@1_indel")
# #
# #     # not in gene
# #     with pytest.raises(Exception):
# #         assert reference.valid_gene_mutation("F@2000_indel")
# #
# # def test_Genome_valid_genome_variant():
# #
# #     assert reference.valid_genome_variant("1a>c")
# #     assert reference.valid_genome_variant("2c>g")
# #     assert reference.valid_genome_variant("3g>a")
# #     assert reference.valid_genome_variant("13335t>g")
# #
# #     # incorrect reference base
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("1t>c")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("2a>g")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("3t>a")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("13335c>g")
# #
# #     # badly formed reference base
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("11c")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("2?>g")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("3P>a")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("13335 >g")
# #
# #     # out of range index
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("0a>c")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("-1c>g")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("-2g>a")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("13336t>g")
# #
# #     # badly formed index
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("1.1a>c")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("2c>fg")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("yg>a")
# #     with pytest.raises(Exception):
# #         assert reference.valid_genome_variant("tg")
# #
# # def test_Genome_convert_variant_to_mutation():
# #
# #     assert reference.convert_variant_to_mutation("N@a1g")=='N@M1V'
# #     assert reference.convert_variant_to_mutation("N@a1c")=='N@M1L'
# #     assert reference.convert_variant_to_mutation("N@a1t")=='N@M1L'
# #     assert reference.convert_variant_to_mutation("N@t2c")=='N@M1T'
# #     assert reference.convert_variant_to_mutation("N@t2a")=='N@M1K'
# #     assert reference.convert_variant_to_mutation("N@t2g")=='N@M1R'
# #     assert reference.convert_variant_to_mutation("N@t2o")=='N@M1O'
# #     assert reference.convert_variant_to_mutation("N@t2x")=='N@M1X'
# #
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@a1g")=='N@M1L'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N1@a1g")=='N@M1V'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@a1g_3")=='N@M1L'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@t1g")=='N@M1L'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@y1g")=='N@M1L'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@a1o")=='N@M1L'
# #     with pytest.raises(Exception):
# #         assert reference.convert_variant_to_mutation("N@a-1g")=='N@M1L'
# #
# # def test_Genome_gbk_fasta_identical():
# #
# #     assert reference.genome_length==reference2.genome_length
# #
# #     assert numpy.array_equal(reference.genome_coding_strand,reference2.genome_coding_strand)
# #
# #     assert numpy.array_equal(reference.genome_index,reference2.genome_index)
# #
# #
# # def test_Genome___repr__():
# #
# #     assert reference.__repr__()=='NC_004148.2\nHuman metapneumovirus\nHMPV\n13335 bases\nacg...cgt'
# #
# #
# # def test_Genome___sub__():
# #
# #     sample=copy.deepcopy(reference)
# #
# #     sample.genome_coding_strand[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based
# #
# #     indices=reference-sample
# #
# #     assert indices[0]==3
# #
# #     assert sample.genome_coding_strand[2]=='t'
# #
# # def test_Genome_contains_gene():
# #
# #     assert reference.contains_gene("F")==True
# #     assert reference.contains_gene("FG")==False
# #     assert reference.contains_gene(5)==False
# #     assert reference.contains_gene("")==False
# #     assert reference.contains_gene("rpoBC")==False
# #     assert reference.contains_gene("RPOB")==False
# #
# # def test_Genome_at_index():
# #
# #     assert reference.at_index(4686)=='F'
# #     assert reference.at_index(4687)=='M2'
# #     assert reference.at_index(5293)=='M2_2'
# #     assert reference.at_index(5450)=='M2_2'
# #     assert reference.at_index(5451)=='SH'
# #     assert reference.at_index(7032) is None
# #     assert reference.at_index(7033)=="L"
# #     assert reference.at_index(13150)=="L"
# #     assert reference.at_index(13151) is None
# #
# #     # wrong gene
# #     assert reference.at_index(7033)!="F"
# #     assert reference.at_index(7033) is not None
# #
# #     # bad position
# #     with pytest.raises(Exception):
# #         reference.at_index(-2)
# #     with pytest.raises(Exception):
# #         reference.at_index(0)
# #     with pytest.raises(Exception):
# #         reference.at_index(1.3)
# #     with pytest.raises(Exception):
# #         reference.at_index('gh')
# #     with pytest.raises(Exception):
# #         reference.at_index(13336)
# #
# # def test_Genome_calculate_snp_distance():
# #
# #     sample=copy.deepcopy(reference)
# #
# #     sample.genome_coding_strand[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based
# #     assert sample.snp_distance(reference)==1
# #
# #     # reverse the change
# #     sample.genome_coding_strand[2]='g'
# #     assert sample.snp_distance(reference)==0
# #
# #     # now change two bases
# #     sample.genome_coding_strand[2]='t'
# #     sample.genome_coding_strand[3]='t'
# #     assert sample.snp_distance(reference)==2
# #
# # sample_01=copy.deepcopy(reference)
# # sample_01.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'],total_coverage_threshold=5,metadata_thresholds={'GT_CONF':5})
# #
# # def test_Genome_apply_vcf():
# #
# #     indices=reference-sample_01
# #     assert indices[0]==4687
# #     assert reference.genome_coding_strand[reference.genome_index==indices[0]]=='t'
# #     assert sample_01.genome_coding_strand[reference.genome_index==indices[0]]=='c'
# #     assert sample_01.coverage[sample_01.genome_index==4687]==68
# #     assert sample_01.genome_sequence_metadata['GT_CONF'][sample_01.genome_index==4687][0]==pytest.approx(613.77)
# #
# #
# #     assert indices[1]==4725
# #     assert reference.genome_coding_strand[reference.genome_index==indices[1]]=='t'
# #     assert sample_01.genome_coding_strand[reference.genome_index==indices[1]]=='c'
# #     assert sample_01.coverage[sample_01.genome_index==4725]==68
# #     assert sample_01.genome_sequence_metadata['GT_CONF'][sample_01.genome_index==4725][0]==pytest.approx(613.77)
# #
# # def test_Genome_list_variants_wrt():
# #
# #     assert sample_01.list_variants_wrt(reference)==['4687t>c','4725t>c','4730c>z','13333c>z','4735_indel','4740_indel']
# #
# #
# # def test_Genome_table_variants_wrt():
# #
# #     # assert sample_01.table_variants_wrt(reference)==['4687t>c']
# #     foo=sample_01.table_variants_wrt(reference)
# #     pass
# #
# #
# # def test_Genome__complement():
# #
# #     test_sequence=numpy.array(['a','c','t','g','z','x','o'])
# #
# #     assert numpy.array_equal(reference._complement(test_sequence),numpy.array(['t','g','a','c','z','x','o']))
# #
# # sample_02=copy.deepcopy(reference)
# # sample_02.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=False,metadata_fields=['GT_CONF'],focussed_indices={4724,4725,4726})
# #
# # def test_sample_02_list_variants_wrt():
# #
# #     assert sample_02.list_variants_wrt(reference)==['4687t>c','4725t>o','4730c>z','4735_indel']
# #
# # sample_03=copy.deepcopy(reference)
# # sample_03.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=False,metadata_fields=['GT_CONF'],focussed_indices={4724,4725,4726,13149,13333})
# #
# # def test_sample_03_list_variants_wrt():
# #
# #     assert sample_03.list_variants_wrt(reference)==['4687t>c','4725t>o','4730c>z','13149g>x','13333c>o','4735_indel']
# #
# #
# # sample_04=copy.deepcopy(reference)
# # sample_04.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'],focussed_indices={4724,4725,4726,13149,13333},mask_file=TEST_CASE_DIR+"01.all.bed")
# #
# # def test_sample_04_list_variants_wrt():
# #
# #     assert sample_04.list_variants_wrt(reference)==[]
# #
# # sample_05=copy.deepcopy(reference)
# # sample_05.apply_vcf_file(vcf_file=TEST_CASE_DIR+"01.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'],focussed_indices={4724,4725,4726,13149,13333},mask_file=TEST_CASE_DIR+"01.M2.bed")
# #
# # def test_sample_05_list_variants_wrt():
# #
# #     assert sample_05.list_variants_wrt(reference)==['13149g>x','13333c>z','13335t>a']
# #
# # # use the subset argument to speed up the genome creation
# # h37rv=Genome(genbank_file="config/H37rV_v3.gbk",name="H37rV_v3",gene_subset=['katG','rpoB','gyrA'])
# #
# # def test_Genome_apply_fasta():
# #
# #     # use gene M2 which spans 4687-5234 incl.
# #     # if we use sample_01 there are variants 4687t>c, 4725t>c already in gene M2 and indels 4735_ins_2 and 4740_ins_2
# #     # hence try applying a fasta AFTER the vcf which contains Ns that starts after the first M2 indel (which is after both SNPs) and finishs in the following M2_2 gene
# #
# #     sample_01.apply_fasta(fasta_file=TEST_CASE_DIR+'01.fasta')
# #
# #     # check that the indel has been flagged
# #     assert sample.is_indel[sample.genome_index==5000]
# #
# #     # check that the indel length has been recorded (note it is negative since the Ns are a deletion)
# #     assert sample.indel_length[sample.genome_index==5000]==-300
# #
# #     # the indel should then be picked up and identified
# #     assert sample_01.list_variants_wrt(reference)==['4687t>c','4725t>c','4730c>z','13333c>z','4735_indel','4740_indel','5000_indel']
# #
# # def test_Genome_H37rV_katG():
# #
# #     reversed_complement=' '.join(h37rv.genome_sequence[h37rv.genome_feature_name=='katG'])[::-1]
# #
# #     # below string taken from mycobrowser.epfl.ch and includes 100 bases upstream
# #     assert reversed_complement =='gtcatctactggggtctatgtcctgattgttcgatatccgacacttcgcgatcacatccgtgatcacagcccgataacaccaactcctggaaggaatgctgtgcccgagcaacacccacccattacagaaaccaccaccggagccgctagcaacggctgtcccgtcgtgggtcatatgaaataccccgtcgagggcggcggaaaccaggactggtggcccaaccggctcaatctgaaggtactgcaccaaaacccggccgtcgctgacccgatgggtgcggcgttcgactatgccgcggaggtcgcgaccatcgacgttgacgccctgacgcgggacatcgaggaagtgatgaccacctcgcagccgtggtggcccgccgactacggccactacgggccgctgtttatccggatggcgtggcacgctgccggcacctaccgcatccacgacggccgcggcggcgccgggggcggcatgcagcggttcgcgccgcttaacagctggcccgacaacgccagcttggacaaggcgcgccggctgctgtggccggtcaagaagaagtacggcaagaagctctcatgggcggacctgattgttttcgccggcaactgcgcgctggaatcgatgggcttcaagacgttcgggttcggcttcggccgggtcgaccagtgggagcccgatgaggtctattggggcaaggaagccacctggctcggcgatgagcgttacagcggtaagcgggatctggagaacccgctggccgcggtgcagatggggctgatctacgtgaacccggaggggccgaacggcaacccggaccccatggccgcggcggtcgacattcgcgagacgtttcggcgcatggccatgaacgacgtcgaaacagcggcgctgatcgtcggcggtcacactttcggtaagacccatggcgccggcccggccgatctggtcggccccgaacccgaggctgctccgctggagcagatgggcttgggctggaagagctcgtatggcaccggaaccggtaaggacgcgatcaccagcggcatcgaggtcgtatggacgaacaccccgacgaaatgggacaacagtttcctcgagatcctgtacggctacgagtgggagctgacgaagagccctgctggcgcttggcaatacaccgccaaggacggcgccggtgccggcaccatcccggacccgttcggcgggccagggcgctccccgacgatgctggccactgacctctcgctgcgggtggatccgatctatgagcggatcacgcgtcgctggctggaacaccccgaggaattggccgacgagttcgccaaggcctggtacaagctgatccaccgagacatgggtcccgttgcgagataccttgggccgctggtccccaagcagaccctgctgtggcaggatccggtccctgcggtcagccacgacctcgtcggcgaagccgagattgccagccttaagagccagatccgggcatcgggattgactgtctcacagctagtttcgaccgcatgggcggcggcgtcgtcgttccgtggtagcgacaagcgcggcggcgccaacggtggtcgcatccgcctgcagccacaagtcgggtgggaggtcaacgaccccgacggggatctgcgcaaggtcattcgcaccctggaagagatccaggagtcattcaactccgcggcgccggggaacatcaaagtgtccttcgccgacctcgtcgtgctcggtggctgtgccgccatagagaaagcagcaaaggcggctggccacaacatcacggtgcccttcaccccgggccgcacggatgcgtcgcaggaacaaaccgacgtggaatcctttgccgtgctggagcccaaggcagatggcttccgaaactacctcggaaagggcaacccgttgccggccgagtacatgctgctcgacaaggcgaacctgcttacgctcagtgcccctgagatgacggtgctggtaggtggcctgcgcgtcctcggcgcaaactacaagcgcttaccgctgggcgtgttcaccgaggcctccgagtcactgaccaacgacttcttcgtgaacctgctcgacatgggtatcacctgggagccctcgccagcagatgacgggacctaccagggcaaggatggcagtggcaaggtgaagtggaccggcagccgcgtggacctggtcttcgggtccaactcggagttgcgggcgcttgtcgaggtctatggcgccgatgacgcgcagccgaagttcgtgcaggacttcgtcgctgcctgggacaaggtgatgaacctcgacaggttcgacgtgcgctga'
# #
# #     indices=h37rv.genome_index[h37rv.genome_feature_name=='katG']
# #     assert indices[0]==2153889
# #     assert indices[-1]==2156211
# #
# #     positions=h37rv.genome_nucleotide_number[h37rv.genome_feature_name=='katG']
# #
# #     # check there is no zero position
# #     assert numpy.sum(positions==0)==0
# #     assert positions[0]==2223
# #     assert positions[-1]==-100
# #
# #     # aminoacids=' '.join(i for i in )
# #
# # def test_Genome_H37rV_rpoB():
# #
# #     sequence=' '.join(h37rv.genome_sequence[h37rv.genome_feature_name=='rpoB'])
# #
# #     # below string taken from mycobrowser.epfl.ch and includes 100 bases upstream
# #     assert sequence =='cgccggccgaaaccgacaaaattatcgcggcgaacgggcccgtgggcaccgctcctctaagggctctcgttggtcgcatgaagtgctggaaggatgcatcttggcagattcccgccagagcaaaacagccgctagtcctagtccgagtcgcccgcaaagttcctcgaataactccgtacccggagcgccaaaccgggtctccttcgctaagctgcgcgaaccacttgaggttccgggactccttgacgtccagaccgattcgttcgagtggctgatcggttcgccgcgctggcgcgaatccgccgccgagcggggtgatgtcaacccagtgggtggcctggaagaggtgctctacgagctgtctccgatcgaggacttctccgggtcgatgtcgttgtcgttctctgaccctcgtttcgacgatgtcaaggcacccgtcgacgagtgcaaagacaaggacatgacgtacgcggctccactgttcgtcaccgccgagttcatcaacaacaacaccggtgagatcaagagtcagacggtgttcatgggtgacttcccgatgatgaccgagaagggcacgttcatcatcaacgggaccgagcgtgtggtggtcagccagctggtgcggtcgcccggggtgtacttcgacgagaccattgacaagtccaccgacaagacgctgcacagcgtcaaggtgatcccgagccgcggcgcgtggctcgagtttgacgtcgacaagcgcgacaccgtcggcgtgcgcatcgaccgcaaacgccggcaaccggtcaccgtgctgctcaaggcgctgggctggaccagcgagcagattgtcgagcggttcgggttctccgagatcatgcgatcgacgctggagaaggacaacaccgtcggcaccgacgaggcgctgttggacatctaccgcaagctgcgtccgggcgagcccccgaccaaagagtcagcgcagacgctgttggaaaacttgttcttcaaggagaagcgctacgacctggcccgcgtcggtcgctataaggtcaacaagaagctcgggctgcatgtcggcgagcccatcacgtcgtcgacgctgaccgaagaagacgtcgtggccaccatcgaatatctggtccgcttgcacgagggtcagaccacgatgaccgttccgggcggcgtcgaggtgccggtggaaaccgacgacatcgaccacttcggcaaccgccgcctgcgtacggtcggcgagctgatccaaaaccagatccgggtcggcatgtcgcggatggagcgggtggtccgggagcggatgaccacccaggacgtggaggcgatcacaccgcagacgttgatcaacatccggccggtggtcgccgcgatcaaggagttcttcggcaccagccagctgagccaattcatggaccagaacaacccgctgtcggggttgacccacaagcgccgactgtcggcgctggggcccggcggtctgtcacgtgagcgtgccgggctggaggtccgcgacgtgcacccgtcgcactacggccggatgtgcccgatcgaaacccctgaggggcccaacatcggtctgatcggctcgctgtcggtgtacgcgcgggtcaacccgttcgggttcatcgaaacgccgtaccgcaaggtggtcgacggcgtggttagcgacgagatcgtgtacctgaccgccgacgaggaggaccgccacgtggtggcacaggccaattcgccgatcgatgcggacggtcgcttcgtcgagccgcgcgtgctggtccgccgcaaggcgggcgaggtggagtacgtgccctcgtctgaggtggactacatggacgtctcgccccgccagatggtgtcggtggccaccgcgatgattcccttcctggagcacgacgacgccaaccgtgccctcatgggggcaaacatgcagcgccaggcggtgccgctggtccgtagcgaggccccgctggtgggcaccgggatggagctgcgcgcggcgatcgacgccggcgacgtcgtcgtcgccgaagaaagcggcgtcatcgaggaggtgtcggccgactacatcactgtgatgcacgacaacggcacccggcgtacctaccggatgcgcaagtttgcccggtccaaccacggcacttgcgccaaccagtgccccatcgtggacgcgggcgaccgagtcgaggccggtcaggtgatcgccgacggtccctgtactgacgacggcgagatggcgctgggcaagaacctgctggtggccatcatgccgtgggagggccacaactacgaggacgcgatcatcctgtccaaccgcctggtcgaagaggacgtgctcacctcgatccacatcgaggagcatgagatcgatgctcgcgacaccaagctgggtgcggaggagatcacccgcgacatcccgaacatctccgacgaggtgctcgccgacctggatgagcggggcatcgtgcgcatcggtgccgaggttcgcgacggggacatcctggtcggcaaggtcaccccgaagggtgagaccgagctgacgccggaggagcggctgctgcgtgccatcttcggtgagaaggcccgcgaggtgcgcgacacttcgctgaaggtgccgcacggcgaatccggcaaggtgatcggcattcgggtgttttcccgcgaggacgaggacgagttgccggccggtgtcaacgagctggtgcgtgtgtatgtggctcagaaacgcaagatctccgacggtgacaagctggccggccggcacggcaacaagggcgtgatcggcaagatcctgccggttgaggacatgccgttccttgccgacggcaccccggtggacattattttgaacacccacggcgtgccgcgacggatgaacatcggccagattttggagacccacctgggttggtgtgcccacagcggctggaaggtcgacgccgccaagggggttccggactgggccgccaggctgcccgacgaactgctcgaggcgcagccgaacgccattgtgtcgacgccggtgttcgacggcgcccaggaggccgagctgcagggcctgttgtcgtgcacgctgcccaaccgcgacggtgacgtgctggtcgacgccgacggcaaggccatgctcttcgacgggcgcagcggcgagccgttcccgtacccggtcacggttggctacatgtacatcatgaagctgcaccacctggtggacgacaagatccacgcccgctccaccgggccgtactcgatgatcacccagcagccgctgggcggtaaggcgcagttcggtggccagcggttcggggagatggagtgctgggccatgcaggcctacggtgctgcctacaccctgcaggagctgttgaccatcaagtccgatgacaccgtcggccgcgtcaaggtgtacgaggcgatcgtcaagggtgagaacatcccggagccgggcatccccgagtcgttcaaggtgctgctcaaagaactgcagtcgctgtgcctcaacgtcgaggtgctatcgagtgacggtgcggcgatcgaactgcgcgaaggtgaggacgaggacctggagcgggccgcggccaacctgggaatcaatctgtcccgcaacgaatccgcaagtgtcgaggatcttgcgtaa'
# #
# #     indices=h37rv.genome_index[h37rv.genome_feature_name=='rpoB']
# #     assert indices[0]==759707
# #     assert indices[-1]==763325
# #
# #     positions=h37rv.genome_nucleotide_number[h37rv.genome_feature_name=='rpoB']
# #
# #     # check there is no zero position
# #     assert numpy.sum(positions==0)==0
# #     assert positions[0]==-100
# #     assert positions[-1]==3519
# #
# #
# #     # sample_04=copy.deepcopy(reference)
# #     # sample_04.apply_vcf_file(vcf_file=TEST_CASE_DIR+"04.vcf",ignore_status=True,ignore_filter=True,metadata_fields=['GT_CONF'])
# #     # (original_bases,indices,new_bases)=reference-sample_04
# #     # assert original_bases[0]=='c'
# #     # assert indices[0]==2155168
# #     # assert new_bases[0]=='g'
# #     # assert sample_04.coverage[sample_04.index==2155168]==53
# #     # assert sample_04.sequence_metadata['GT_CONF'][sample_04.index==2155168][0]==pytest.approx(500.23)
# #
# #
# # def test_Genome_H37rV_site_05():
# #
# #     sample_tb_01=copy.deepcopy(h37rv)
# #     sample_tb_01.apply_vcf_file(vcf_file="tests/test-cases/05.vcf",\
# #                                ignore_status=True,ignore_filter=False,\
# #                                metadata_fields=['GT_CONF','GT_CONF_PERCENTILE','FRS','DPF','DP'])
# #
# #     assert sample_tb_01.list_variants_wrt(h37rv)==['7362g>c', '9304g>a', '761110a>t', '763031t>c', '2154724g>t', '2155168g>c']
# #
# #     assert sample_tb_01.genes['gyrA'].list_mutations_wrt(h37rv.genes['gyrA'])==['E21Q', 'G668D']
# #     assert sample_tb_01.genes['rpoB'].list_mutations_wrt(h37rv.genes['rpoB'])==['D435V', 'A1075A']
# #     assert sample_tb_01.genes['katG'].list_mutations_wrt(h37rv.genes['katG'])==['S315T', 'R463L']
# #
# #
# # def test_Genome_H37rV_site_05_mask():
# #
# #     # the mask should exclude the katG@S315T mutation
# #
# #     sample_tb_01=copy.deepcopy(h37rv)
# #     sample_tb_01.apply_vcf_file(vcf_file="tests/test-cases/05.vcf",\
# #                                ignore_status=True,ignore_filter=False,\
# #                                metadata_fields=['GT_CONF','GT_CONF_PERCENTILE','FRS','DPF','DP'],\
# #                                mask_file="tests/test-cases/05.bed")
# #
# #     assert sample_tb_01.list_variants_wrt(h37rv)==['7362g>c', '9304g>a', '761110a>t', '763031t>c', '2154724g>t']
# #
# #     assert sample_tb_01.genes['gyrA'].list_mutations_wrt(h37rv.genes['gyrA'])==['E21Q', 'G668D']
# #     assert sample_tb_01.genes['rpoB'].list_mutations_wrt(h37rv.genes['rpoB'])==['D435V', 'A1075A']
# #     assert sample_tb_01.genes['katG'].list_mutations_wrt(h37rv.genes['katG'])==['R463L']
# #
# #
# # def test_Genome_H37rV_site_05_focussed_indices():
# #
# #     # be supplying a set of genome indices the code will return null and filter fail calls at those positions
# #     mask=numpy.isin(h37rv.genome_feature_name,['gyrA'])
# #     focussed_indices=set(h37rv.genome_index[mask])
# #
# #     sample_tb_01=copy.deepcopy(h37rv)
# #     sample_tb_01.apply_vcf_file(vcf_file="tests/test-cases/05.vcf",\
# #                                ignore_status=True,ignore_filter=False,\
# #                                metadata_fields=['GT_CONF','GT_CONF_PERCENTILE','FRS','DPF','DP'],\
# #                                focussed_indices=focussed_indices)
# #
# #     assert sample_tb_01.list_variants_wrt(h37rv)==['7362g>c','7547c>o','7550c>x','7553c>x','7556c>x','7559g>x','7562c>x','9304g>a','761110a>t','763031t>c','2154724g>t','2155168g>c']
# #
# #     assert sample_tb_01.genes['gyrA'].list_mutations_wrt(h37rv.genes['gyrA'])==['E21Q', 'G82O', 'N83X', 'Y84X', 'H85X', 'P86X', 'H87X', 'G668D']
# #     assert sample_tb_01.genes['rpoB'].list_mutations_wrt(h37rv.genes['rpoB'])==['D435V', 'A1075A']
# #     assert sample_tb_01.genes['katG'].list_mutations_wrt(h37rv.genes['katG'])==['S315T', 'R463L']
