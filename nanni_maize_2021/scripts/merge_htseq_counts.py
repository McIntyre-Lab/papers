#!/usr/bin/env python3

# Take in 2 results files from merged sample htseq-count PE and SE alignments and sum them

import pandas as pd
import argparse
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Take in 2 sample merged results files from htseq-count PE and SE alignments and sum them")

    # Input arguments
    parser.add_argument("-p", "--paired", dest="inPE", required=True, help="Results file from PE htseq-count output after merging all samples (1 per column)")
    parser.add_argument("-s", "--single", dest="inSE", required=True, help="Results file from SE htseq-count output after merging all samples (1 per column)")

    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output filename for the summed expression counts")

    args = parser.parse_args()
    return args

def main():
    # Check that input files exist
    if not os.path.isfile(args.inSE) or not os.path.isfile(args.inPE):
        raise FileNotFoundError

    # Open input files
    PE = pd.read_csv(args.inPE,sep="\t")
    SE = pd.read_csv(args.inSE,sep="\t")

    # Stack files and group by transcript_id, sum counts
    total = PE.append(SE).groupby('gene_id').sum().reset_index()

    # Output summed file
    total.to_csv(args.outFile,sep="\t",index=False)


if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()
