#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get frequencies of fragment sex bias flags")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of fragment-level annotation flags")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file for frequencies of sex bias gene expression characterizations")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    annotDF = pd.read_csv(args.inFile, low_memory=False)

    # Select fragment features only
    fragDF = annotDF.loc[annotDF['featureType']=="fragment"]
    del(annotDF)
    
    # Get frequencies of fragment flags of:
    #   1) ttest significance flags with ratio != 1
    #   2) fold change >= 2
    #   3) ttest significant flag and fold change >= 2
    of = open(args.outFile, 'w')
    of.write("Get frequencies of fragment flags of:\n1) ttest significance flags with ratio != 1\n{}\n".format(
            pd.crosstab(fragDF['ratio_trend'],fragDF['flag_ttest_pval'], dropna=False)))
    of.write("\n2) fold change >= 2\n{}\n".format(
            fragDF['ratio_expressed'].value_counts().to_string()))
    of.write("\n3) ttest significant flag and fold change >= 2\n{}\n".format(
            pd.crosstab(fragDF['ratio_expressed'], fragDF['flag_ttest_pval'])))
    of.close()
    
    # Verify the frequencies match the combo flag counts
    if fragDF['flag_ttest_1_ratio2_M'].sum() != len(fragDF[(fragDF['ratio_expressed']=="male")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    if fragDF['flag_ttest_1_ratio2_F'].sum() != len(fragDF[(fragDF['ratio_expressed']=="fem")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    if fragDF['flag_ttest_1_ratio2_U'].sum() != len(fragDF[(fragDF['ratio_expressed']=="unb")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    if fragDF['flag_ttest_1_trend_M'].sum() != len(fragDF[(fragDF['ratio_trend']=="male")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    if fragDF['flag_ttest_1_trend_F'].sum() != len(fragDF[(fragDF['ratio_trend']=="fem")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    if fragDF['flag_ttest_1_trend_U'].sum() != len(fragDF[(fragDF['ratio_trend']=="unb")&(fragDF['flag_ttest_pval']==1)]) :
        print("WARNING: Flag value frequencies not matching")
    # All match in both mel and sim
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()