#! /usr/bin/python3
"""Load a gumpy Genome object, and use pickle to save to disk.
Due to the security implications of the pickle module, USE AT YOUR OWN RISK. 
It is recommended to avoid sending and recieving pickled objects.
"""
import argparse
import pickle
import os

import gumpy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genbank_file", required=True, help="Path to a single genbank file"
    )
    parser.add_argument("--path", required=True, help="Path to the output file")
    parser.add_argument(
        "--overwrite",
        required=False,
        default=False,
        action="store_true",
        help="Overwrite a file at the given path (if exists)",
    )
    options = parser.parse_args()

    # instantiate a gumpy Genome object using the provided Genbank file
    reference_genome = gumpy.Genome(options.genbank_file, show_progress_bar=True)

    if os.path.exists(options.path):
        if options.overwrite:
            print("File existed! Overwriting")
            pickle.dump(reference_genome, open(options.path, "wb"))
        else:
            print("File already exists! Use --overwrite to replace it")
    else:
        pickle.dump(reference_genome, open(options.path, "wb"))
