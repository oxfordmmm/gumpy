#!/usr/bin/python3
"""Create a catalogue compatable with piezo from a given VCF.
Assumes that all mutations are in the `S` category out of simplicity.
Takes arguments from the command line

Args:
    genome_path (str): Path to the reference genome
    vcf_path (str): Path to the VCF file
    output_filename (str): Desired path to the output csv"""
import gumpy
import sys
import functools


def mutations_to_catalogue(mutations, filename, ref_genome, catalogue):
    """Convert a list of mutations to a catalogue which can be interpreted by piezo

    Args:
        mutations (list(str)): List of mutations in GARC
        filename (str): Path to the output file
        ref_genome (str): Name of the genome
        catalogue (str): Name of the catalogue
    """
    header = (
        "GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,"
        "PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n"
    )
    common_all = f"{ref_genome},{catalogue},1.0,GARC1,RUS,NAN,"

    # Get all genes within the mutations
    genes = sorted(list(set([mutation.split("@")[0] for mutation in mutations])))
    with open(filename, "w") as f:
        f.write(header)
        # Piezo requires all prediction values are in the catalogue,
        # and all valid `PREDICTION_VALUES` include `R`
        # So add a dummy record with R as S and U are used elsewhere
        f.write(common_all + "NOTAGENE@A1!,R,{},{},{}\n")
        # Add general rules for genes
        for gene in genes:
            f.write(common_all + gene + "@*=,S,{},{},{}\n")
            f.write(common_all + gene + "@*?,U,{},{},{}\n")
        # Add actual mutations
        for mutation in mutations:
            line = common_all + mutation + ",S,{},{},{}\n"
            f.write(line)


if __name__ == "__main__":
    assert len(sys.argv) == 4, (
        "Incorrect usage. Try: python to_piezo_catalogue.py "
        "<genome_path> <vcf_path> <output_filename>"
    )
    genome_path = sys.argv[1]
    vcf_path = sys.argv[2]
    output_filename = sys.argv[3]

    # Get the GeneDifference objects to pull out mutations
    g = gumpy.Genome(genome_path)
    testDNA = gumpy.VCFFile(vcf_path)
    diff = testDNA.difference(g)
    gene_diffs = diff.gene_differences()
    # Get a single list of mutations by concatenating the mutations lists
    mutations = functools.reduce(
        lambda x, y: x + y, [g.mutations.tolist() for g in gene_diffs], []
    )
    # Write the catalogue
    mutations_to_catalogue(
        mutations, output_filename, g.name, vcf_path.split("/")[-1].replace(".vcf", "")
    )
