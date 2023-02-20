#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add flags/variables for lists of genes to ChIP RNA gene-level and ortholog files")

    # Input data
    parser.add_argument("-t", "--mel-ttest", dest="inTTEST", required=True, help="Input CSV of D. melanogaster gene-level file with combo flag gene_sex_bias_ttest_foldchange")
    parser.add_argument("-ot", "--ortho-ttest", dest="inOT", required=True, help="Input CSV of D. melangoster and D. simulans orthologous gene-level file with ttest+foldchange values")
    parser.add_argument("-ml", "--mel-gene-list", dest="inMlist", required=True, help="Input CSV of D. melanogaster gene list flags")

    # Output data
    parser.add_argument("-s", "--output-suffix", dest="outSuff", required=True, help="Suffix to use for output files - do not include file extension or underbar before/after suffix (will add to the end of the input CSV file names)")

    args = parser.parse_args()
    return args

def main():
    # Get inputs
    ttestDF = pd.read_csv(args.inTTEST, low_memory=False)
    otDF = pd.read_csv(args.inOT, low_memory=False)
    melList = pd.read_csv(args.inMlist, low_memory=False)
    
    # Rename ttest gene flag to match other files
    ttestDF = ttestDF.rename(columns={'fbgn':'FBgn'})

    # Merge in flags of previous literature gene lists to data files
    mergeTtestDF = pd.merge(
            ttestDF,
            melList,
            how="outer",
            left_on="FBgn",
            right_on="gene_id",
            indicator="merge_check",
            validate="1:1"
    )
    mergeTtestDF = mergeTtestDF[mergeTtestDF["merge_check"]!="right_only"]
    mergeTtestDF = mergeTtestDF.drop(columns=["merge_check"])

    mergeOTDF = pd.merge(
            otDF,
            melList,
            how="outer",
            left_on="mel_geneID",
            right_on="gene_id",
            indicator="merge_check",
            validate="m:1"
    )
    mergeOTDF = mergeOTDF[mergeOTDF["merge_check"]!="right_only"]
    mergeOTDF = mergeOTDF.drop(columns=["merge_check"]) 

    ## Output final merged files
    mergeTtestDF.to_csv(args.inTTEST.split(".csv")[0]+"_"+args.outSuff+".csv", index=False)
    mergeOTDF.to_csv(args.inOT.split(".csv")[0]+"_"+args.outSuff+".csv", index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
