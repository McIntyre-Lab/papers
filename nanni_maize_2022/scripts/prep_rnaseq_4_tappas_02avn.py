#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Prepare files for tappAS (genotype expression matrix of expected counts and design files) for maize ozone RSEM")

    # Input data
    parser.add_argument("-t", "--input-TPM", dest="inTPM", required=True, help="Input TSV of combined rsem expression matrix with TPM on/off flags (used to determine detected transcripts)")
    parser.add_argument("-c", "--input-count", dest="inCount", required=True, help="Input TSV of combined rsem expression matrix with expected counts (used to output expected counts for detected transcripts for tappas)")
    parser.add_argument("-e", "--exclude", dest="exclude", required=False, action='append', help="Samples to exclude from expression matrices, multiple values can be listed with each '-e'")
    parser.add_argument("-v", "--value-type", dest="inVal", required=False, choices=['TPM','expected_count'], default='TPM', help="Value type to include in tappas expression matrix (TPM or expected_count, default: TPM)")

    # Output data
    parser.add_argument("-o", "--output-directory", dest="outDir", required=True, help="Output directory")

    args = parser.parse_args()
    return args

def main():
    # Get combined expression matrix with on/off flags
    tpmDF = pd.read_csv(args.inTPM, sep="\t")
    countDF = pd.read_csv(args.inCount, sep="\t")
    
    # Remove excluded columns if provided
    if args.exclude is not None:
        for s in args.exclude:
            tpmDF = tpmDF.drop(columns=[c for c in tpmDF.columns if c==s])   
            countDF = countDF.drop(columns=[c for c in countDF.columns if c==s])   
            print("Removed {} columns from matrix...".format(s))

    # Count and drop transcripts that are not detected in any samples
    #     and transcripts only detected in one samples
    print("{} transcripts detected in 0 samples\n{} transcripts detected in 1 sample\n...Dropped from tappas files".format(
            len(tpmDF[tpmDF['sum_flag']==0]),len(tpmDF[tpmDF['sum_flag']==1])))
    detectDF = tpmDF[tpmDF['sum_flag']>1].copy()
    
    # Merge detected transcripts wtih expected count values
    # Select only those in both (no left_only so it is everyhting in detected DF)
    if args.inVal == "expected_count":
        mergeDF = pd.merge(detectDF['transcript_id'],countDF,how='outer',on='transcript_id',indicator='merge_check')
        mergeDF = mergeDF[mergeDF['merge_check']=="both"]
    else:
        mergeDF = detectDF.copy()
    
    # Get Amb vs. Oz expression matrices and design files
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        cols = [c for c in mergeDF.columns if (genotype in c) and ('flag' not in c) and ('mean' not in c)]
        mergeDF[['transcript_id']+cols].rename(columns={'transcript_id':""}).to_csv(
                "{}/sbys_{}_4_tappas_{}.tsv".format(args.outDir,genotype,args.inVal),sep="\t",index=False)
        designDF = pd.DataFrame({'Sample':cols})
        designDF['Condition'] = np.where(designDF['Sample'].str.contains('Amb'),"Ambient",
                np.where(designDF['Sample'].str.contains('Ele'),"Ozone","oops"))
        if len(designDF[designDF['Condition']=="oops"])>0:
            print("ERROR: Cannot assign Condition in {}...".format(genotype))
            sys.exit()
        else:
            designDF.sort_values('Condition').to_csv("{}/df_{}_4_tappas.tsv".format(args.outDir,genotype),sep="\t",index=False)
    # Get B73 Amb vs. all  other genotype Amb expression matrices and design files
    B73ambCols = [c for c in mergeDF.columns if ('B73' in c) and ('Amb' in c) and ('flag' not in c) and ('mean' not in c)]
    for genotype in ["C123","Hp301","Mo17","NC338"]:
        AmbCols = [c for c in mergeDF.columns if (genotype in c) and ('Amb' in c) and ('flag' not in c) and ('mean' not in c)]
        mergeDF[['transcript_id']+B73ambCols+AmbCols].rename(columns={'transcript_id':""}).to_csv(
            "{}/sbys_B73_vs_{}_Amb_4_tappas_{}.tsv".format(args.outDir,genotype,args.inVal),sep="\t",index=False)
        designDF = pd.DataFrame({'Sample':B73ambCols+AmbCols})
        designDF['Condition'] = np.where(designDF['Sample'].str.contains('B73'),"B73",
                np.where(designDF['Sample'].str.contains(genotype),genotype,"oops"))
        if len(designDF[designDF['Condition']=="oops"])>0:
            print("ERROR: Cannot assign Condition in {}...".format(genotype))
            sys.exit()
        else:
            designDF.sort_values('Condition').to_csv("{}/df_B73_vs_{}_Amb_4_tappas.tsv".format(args.outDir,genotype),sep="\t",index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

