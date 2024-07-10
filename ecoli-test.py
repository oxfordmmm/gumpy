import argparse
import gumpy


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", help="Path to a reference genome")
    parser.add_argument("--vcf", help="Path to a VCF file")
    args = parser.parse_args()

    print("Building reference genome...")
    reference = gumpy.Genome(args.reference, show_progress_bar=True)
    print()
    print("Parsing VCF...")
    vcf = gumpy.VCFFile(
        args.vcf,
        ignore_filter=True,
        minor_population_indices=set(reference.nucleotide_index.tolist()),
    )
    print()

    if len(vcf.minor_populations) == 0:
        print("No minor alleles found in VCF file!")

    sample = reference + vcf
    if len(sample.minority_populations_GARC(reference=reference)) > 0:
        print("Genome-level minor alleles:")
        for minor_population in sample.minority_populations_GARC(reference=reference):
            print(minor_population)
    else:
        print("No minor alleles found at a genome level!")

