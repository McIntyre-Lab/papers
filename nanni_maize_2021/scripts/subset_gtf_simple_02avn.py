#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import sys
import csv

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Subset Reference GTF based on a list of transcript IDs")

    # Input arguments
    parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="GTF file where transcript_id and gene_id are the first two attributes in the attributes column (order does not matter)")
    parser.add_argument("-t", "--list-type", dest="inType", required=True, choices=['transcript_id','gene_id'], default='transcript_id', help="Type of feature in inclusion/exclusion list provided (transcript_id or gene_id, default is transcript_id)")
    parser.add_argument("-i", "--include-list", dest="inInclude", required=False, help="List of transcript or gene IDs to include in output from input GTF (no header) - cannot be used with exclusion list")
    parser.add_argument("-e", "--exclude-list", dest="inExclude", required=False, help="List of transcript or gene IDs to exclude in output from input GTF (no header) - cannot be used with inclusion list")
  
    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file name for subset GTF")

    args = parser.parse_args()
    return args


def main():
    # Get GTF file and count total lines
    gtf = pd.read_csv(args.inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], dtype=str, sep="\t",low_memory=False)

    # Get list ID variable
    listVar = args.inType
    
    # Check that only one list argument was provided
    if args.inInclude is None and args.inExclude is None:
        print("ERROR: No list provided to subset GTF on")
        sys.exit()
    elif args.inInclude is not None and args.inExclude is None:
        include = True
        idDF = pd.read_csv(args.inInclude, names=[listVar])
    elif args.inInclude is None and args.inExclude is not None:
        include = False
        idDF = pd.read_csv(args.inExclude, names=[listVar])
    else:
        print("ERROR: Cannot use both inclusion and exclusion list - only provide one or the other")
        sys.exit()
    
    # Get attribute labels 1 and 2
    gtf['attribute_name_1'] = gtf['attribute'].str.split(" ").str[0]
    gtf['attribute_name_2'] = gtf['attribute'].str.split(" ").str[2]
    
    # Extract transcript_id and gene_id from attribute column
    gtf['transcript_id'] = np.where(gtf['attribute_name_1']=="transcript_id",
                                       gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                       np.where(gtf['attribute_name_2']=="transcript_id",
                                                gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))
    gtf['gene_id'] = np.where(gtf['attribute_name_1']=="gene_id",
                                       gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                       np.where(gtf['attribute_name_2']=="gene_id",
                                                gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))

    if include:
        # If list is for inclusion
        # Extract GTF elements that have a value contained in list of IDs
        filterDF = gtf[gtf[listVar].isin(idDF[listVar])]
    else:
        # If list is for exclusion
        # Extract GTF elements that do NOT have a value contained in list of IDs
        filterDF = gtf[~gtf[listVar].isin(idDF[listVar])]
        
    # Output filtered file
    filterDF[['chr','source','feature','start','end','score','strand','frame',
             'attribute']].to_csv(args.outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

