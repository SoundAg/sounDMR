
import pandas as pd
import math
import os
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("chunk_size", type=int)
parser.add_argument("output_dir", nargs="?", default="./chunks/")
args = parser.parse_args()


def split_by_chunk(input_file, chunk_size, output_dir = "./chunks/"):
    """
    """
    # Get base name
    base_name = os.path.basename(input_file)
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Read the BED file into a data frame
    column_names = ["chromosome", "start", "end", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    bed_df = pd.read_table(input_file, header = None, names = column_names)
    nucls = bed_df["start"].tail(1).values[0]
    chrmname = bed_df["chromosome"].head(1).values[0]
    # Get number of chunks
    total_chunks = math.ceil(nucls/chunk_size)
    # Initialize variables for chunk creation
    lower_val = 0
    upper_val = chunk_size
    for chunk_n in range(1, total_chunks + 1):
        chunk = bed_df[(bed_df["start"] > lower_val) & (bed_df["start"] <= upper_val)]
        # check if it is blank, if so make a all NA line
        if chunk.shape[0] == 0:
            chunk = pd.DataFrame(columns = column_names)
            chunk.loc[0] = [chrmname, lower_val] + ['NA']*10
            chunk.index.name = None
            print(f'ZERO: {chunk_n}')
        chunkname = "chunk" + str(chunk_n) + "_" + str(lower_val) + "_" + str(upper_val)
        chunkdir = os.path.join(output_dir, chunkname)
        os.makedirs(chunkdir, exist_ok=True)
        chunkdir = os.path.join(chunkdir, base_name)
        chunk.to_csv(chunkdir, sep = "\t", header = False, index = False)
        lower_val = lower_val + chunk_size
        upper_val = upper_val + chunk_size
        #print("Chunk " + str(chunk_n) + " created.")


if __name__ == "__main__":
    split_by_chunk(args.input_file, args.chunk_size, args.output_dir)
    sys.exit(0)

