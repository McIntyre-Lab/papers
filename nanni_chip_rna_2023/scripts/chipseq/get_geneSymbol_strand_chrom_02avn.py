#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import os
import gffutils

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add gene symbol and strand information to list of reference geneIDs")

    # Input data
    parser.add_argument("-i", "--input", dest="inList", required=True, help="Input list of geneIDs with NO HEADER (must be contained within reference annotation)")
    parser.add_argument("-g", "--GTF", dest="inGTF", required=True, help="Reference annotation GTF")
    parser.add_argument("-c", "--chrom", dest="getChrom", action='store_true', help="Add chromosome to list and output X vs. Autosome counts")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args

def main():
    # Get list of gene IDs
    geneDF = pd.read_csv(args.inList,names=['geneID'])

    # Get reference GTF file
    # Create gffutils database file if it does not already exist
    if os.path.isfile(args.inGTF+".db"):
        db = gffutils.FeatureDB(args.inGTF+".db")
    else:
        db = gffutils.create_db(args.inFile,args.inGTF+".db")
    
    # Add geneSymbol and strand columns using GTF
    geneDF['geneSymbol'] = geneDF.apply(lambda x: db[x['geneID']].attributes['gene_symbol'][0], axis=1)
    geneDF['strand'] = geneDF.apply(lambda x: db[x['geneID']].strand, axis=1)

    # Add chromosome if argument provided
    if args.getChrom:
        geneDF['chr'] = geneDF.apply(lambda x: db[x['geneID']].chrom, axis=1)
        geneDF['flag_x'] = np.where((geneDF['chr']=="X") | (geneDF['chr']=="Scf_X"),"1",
                        np.where((geneDF['chr'].str.contains("Scaffold"))|(geneDF['chr'].str.contains("NODE")),"-1","0"))
        geneDF['flag_4'] = np.where((geneDF['chr']=="4") | (geneDF['chr']=="Scf_4"),"1",
                        np.where((geneDF['chr'].str.contains("Scaffold"))|(geneDF['chr'].str.contains("NODE")),"-1","0"))
        print("X\t"+str(len(geneDF[geneDF['flag_x']=="1"])))
        print("4\t"+str(len(geneDF[geneDF['flag_4']=="1"])))
        print("Auto\t"+str(len(geneDF[(geneDF['flag_x']=="0")&(geneDF['flag_4']=="0")])))
        print("Contig\t"+str(len(geneDF[geneDF['flag_x']=="-1"])))
        geneDF = geneDF.drop(columns=['chr','flag_x','flag_4'])

    # Check for any missing values
    if len(geneDF[geneDF['geneSymbol'].isna()]) > 0:
        print("WARNING: "+str(len(geneDF[geneDF['geneSymbol'].isna()]))+" GENE SYMBOL VALUES ARE MISSING")
    if len(geneDF[geneDF['strand'].isna()]) > 0:
        print("WARNING: "+str(len(geneDF[geneDF['strand'].isna()]))+" STRAND VALUES ARE MISSING")
    
    # Output list with new columns
    geneDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

