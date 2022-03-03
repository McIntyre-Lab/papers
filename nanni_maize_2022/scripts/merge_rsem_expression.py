#!/usr/bin/env python3

# Take in 2 results files from rsem-calculate-expression PE and SE alignments and sum them

import pandas as pd
import argparse
import os
import sys


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Take in 2 results files from rsem-calculate-expression PE and	SE alignments and sum them")

    # Input arguments
    parser.add_argument("-p", "--paired", dest="inPE", required=True, help="Results file from PE rsem-calculate-expression output")
    parser.add_argument("-s", "--single", dest="inSE", required=True, help="Results file from PE rsem-calculate-expression output")

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

    # Check that columns are in the same order (should be)
    

    # Stack files and group by transcript_id, sum counts
    total = PE.append(SE).groupby('transcript_id').sum().reset_index()

    # Output summed file
    total.to_csv(args.outFile, index=False, sep='\t')


if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()
