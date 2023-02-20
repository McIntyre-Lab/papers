#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add flag to list of genes based on smaller list")

    # Input data
    parser.add_argument("-m", "--main", dest="inMain", required=True, help="Input CSV file of main list of genes (gene_id columns)")
    parser.add_argument("-i", "--input", dest="inInput", required=True, help="Input list of genes to be flagged")
    parser.add_argument("-f", "--flag-name", dest="inName", required=True, help="Input flag name")

    # Output data
    parser.add_argument("-o", "--outFile", dest="outFile", required=True, help="Output CSV file for main file with new flag added")

    args = parser.parse_args()
    return args

def main():
    # Get input files
    mainDF = pd.read_csv(args.inMain).astype(str)
    inDF = pd.read_csv(args.inInput).astype(str)
    print("\t"+str(len(mainDF))+" in main list\n\t"+str(len(inDF))+" in list to be flagged")
    
    # Make new flag
    mainDF[args.inName] = np.where(mainDF['gene_id'].isin(inDF['gene_id']),"1","0")
    
    # Check that all genes in input list are in the main list
    print("\n\t"+str(len(inDF[~inDF['gene_id'].isin(mainDF['gene_id'])]))+" genes not in main list:")
    print(inDF[~inDF['gene_id'].isin(mainDF['gene_id'])].to_csv(index=False,header=False))

    # Output main with new flag
    mainDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

