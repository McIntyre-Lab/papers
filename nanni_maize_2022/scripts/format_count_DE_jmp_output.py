#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Format GO enrichment significance flags")

    # Input data
    parser.add_argument("--GOterm", dest="inTerm", required=True, help="Input GO term TSV file")
    parser.add_argument("--JMPfile", dest="inJMP", required=True, help="")

    # Output data
    parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Output prefix for final CSV files (full table and signifiant flags)")

    args = parser.parse_args()
    return args

def main():
    # Get input files
    inDF = pd.read_csv(args.inJMP)
    termDF = pd.read_csv(args.inTerm, sep="\t", quoting=3, low_memory=False)
    termDF.columns = termDF.columns.str.strip("\"")
    termDF = termDF.apply(lambda x: x.str.strip("\""))
    termDF = termDF[['go_id','Term','Secondary']].drop_duplicates()

    # Add significance flag
    inDF['flag_go_sig'] = np.where(inDF['Fisher_FDR_p']<=0.01,1,0)
    
    # Fix "go" to be "GO"
    inDF['Category'] = "GO" + inDF['Category'].str[2:]
    
    # Merge in GO terms
    mergeDF1 = pd.merge(termDF,inDF,how='outer',left_on="go_id",right_on="Category",indicator="merge_check1")
    mergeDF2 = pd.merge(termDF,mergeDF1[mergeDF1['merge_check1']=="right_only"][inDF.columns],how='outer',left_on="Secondary",right_on="Category",indicator="merge_check2")
    cols = [c for c in mergeDF1.columns if (c in inDF.columns) or (c == "Term")]
    mergeDF = pd.concat([mergeDF1[mergeDF1['merge_check1']=="both"][cols],mergeDF2[mergeDF2['merge_check2']=="both"][cols]]).drop_duplicates().sort_values(by=['SigVarName','Fisher_FDR_p'],axis=0)
    mergeDF['flag_go_sig'] = mergeDF['flag_go_sig'].fillna(0)

    # Counts
    print("\nSignificant GO per flag:\n")
    print(mergeDF.groupby('SigVarName')['flag_go_sig'].sum().to_string()+"\n\n")
    print("\nSignificant flags per GO:\n")
    temp = mergeDF.groupby('Term')['flag_go_sig'].sum()
    print(temp[temp>0].to_string()+"\n\n")

    # Get significant values for each flag, drop columns without any significant flags
    finalDF = mergeDF.pivot_table(index='Term',columns='SigVarName',values='flag_go_sig').fillna(0).astype(int)
    finalDF = finalDF[finalDF.sum(axis=1)>=1].reset_index()
    finalDF.to_csv(args.outPrefix+"_sig_terms.csv",index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
