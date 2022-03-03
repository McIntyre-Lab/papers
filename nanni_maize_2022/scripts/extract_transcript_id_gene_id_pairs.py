#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Extract all unique pairs of individual transcript_id to gene_id")

    # Input data
    parser.add_argument("-i", "--input", dest="inAnnot", required=True, help="CSV file of event_id to transcript_id to gene_id (*_event2transcript2gene_index.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file")

    args = parser.parse_args()
    return args

def main():
    # Get input file of event_id to transcript_id to gene_id
    annotDF = pd.read_csv(args.inAnnot, low_memory=False)

    # Extract all unique pairs of individual transcript_id to gene_id
    # This is the file everything else will be merged on
    xcrptGene = annotDF[(~annotDF['transcript_id'].str.contains("|",regex=False))&
            (~annotDF['gene_id'].str.contains("|",regex=False))&
            (annotDF['transcript_id']!="Unannotated")][['gene_id','transcript_id']].drop_duplicates()

    xcrptGene.to_csv(args.outFile, index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

