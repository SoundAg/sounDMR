"""
Split by chromosome

"""

import os
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()

def split_by_chromosome(input_file):
    """
    input_file: A string with the name of the input file. It uses a
    full (absolute) path.
    Scaffolds, if present, should be at the end of the files.
    return: None
    """
    # Get dir only
    input_dir = os.path.dirname(input_file)
    base_name = os.path.basename(input_file)
    # Open the input file for reading
    with open(input_file, "r") as fo:
        # Create a list to store the file connections for each chromosome output file
        output_files_conn = {}
        # Read the input file line by line and process each line
        while True:
            line = fo.readline()
            if not line:
                break
            # Split the line by tab to get the chromosome
            chromosome = line.split("\t")[0]
            # skip line if chromosome starts with "scaf"
            if chromosome.lower().startswith("scaf"):
                continue
            # Create the output directory for the chromosome if not already created
            if not os.path.exists(os.path.join(input_dir, "chr_" + chromosome)):
                os.makedirs(os.path.join(input_dir, "chr_" + chromosome))
            # Create the output file for the chromosome if not already opened
            if chromosome not in output_files_conn:
                output_file = os.path.join(input_dir, "chr_" + chromosome, base_name)
                output_files_conn[chromosome] = open(output_file, "w")
                print(f'starting chromosome {chromosome}')
            output_files_conn[chromosome].write(line)
    # close all connections inside output_files_conn
    for chromosome in output_files_conn:
        output_files_conn[chromosome].close()
        print(f'closed file handle for chromosome {chromosome}')

if __name__ == "__main__":
    split_by_chromosome(args.path)
    sys.exit(0)
