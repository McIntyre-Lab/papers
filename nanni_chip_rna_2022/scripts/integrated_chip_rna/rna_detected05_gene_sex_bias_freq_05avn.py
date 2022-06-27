#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get frequencies of detection of genes APN > 5 with gene expression characterizations of sex bias")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of gene-level annotation flags with combo flag gene_sex_bias_ttest_foldchange")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file for frequencies of gene expression characterization and APN > 5 detection flag")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    annotDF = pd.read_csv(args.inFile, low_memory=False)
    
    # To ensure nan values are counted in crosstab, convert to string
    annotDF['gene_ratio2_ttest'] = annotDF['gene_ratio2_ttest'].fillna("nan")
    annotDF['gene_trend_ttest'] = annotDF['gene_trend_ttest'].fillna("nan")
    
    # Crosstab of flag_rna_detected05 and gene_sex_bias
    of = open(args.outFile, 'w')
    trendRatioStack = pd.DataFrame(annotDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['flag_rna_detected05'].value_counts()).rename(columns={'flag_rna_detected05':'count'})
    trendRatioWide = trendRatioStack.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='flag_rna_detected05', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))

    of.write("\n{}\n".format(pd.crosstab(annotDF['gene_sex_bias'], annotDF['flag_rna_detected05'])))
    of.write("\n{}\n".format(pd.crosstab(annotDF['gene_trend_ttest'], annotDF['flag_rna_detected05'])))
    of.write("\n{}\n".format(pd.crosstab(annotDF['gene_ratio2_ttest'], annotDF['flag_rna_detected05'])))
    of.write("\n{}\n".format(trendRatioWide))

    # Split by xsome and by sex-limited/sex-bias
    ## Sex-limited (have to check if any exist in each set since there are so few sex-limited)
    limited = annotDF[(annotDF['gene_rna']=="fem")|(annotDF['gene_rna']=="male")]
    if len(limited[limited['flag_rna_detected05']==1]) == 0:
        detectedXsomeWideLimited = pd.DataFrame()
    else:
        detectedXsomeStackLimited = pd.DataFrame(limited[limited['flag_rna_detected05']==1].groupby(
                ['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
        detectedXsomeWideLimited = detectedXsomeStackLimited.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                     columns='xsome', fill_value=0).droplevel(0,axis=1)
    if len(limited[limited['flag_rna_detected05']==0]) == 0:
        noDetectXsomeWideLimited = pd.DataFrame()
    else:
        noDetectXsomeStackLimited = pd.DataFrame(limited[limited['flag_rna_detected05']==0].groupby(
                ['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
        noDetectXsomeWideLimited = noDetectXsomeStackLimited.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                     columns='xsome', fill_value=0).droplevel(0,axis=1)
    of.write("\n#### Sex-limited genes only (gene_rna == male or fem) ####\n\n")
    if not noDetectXsomeWideLimited.empty:
        of.write("\nGenes not APN>5 (flag_rna_detected05==0) by xsome\n{}\n".format(noDetectXsomeWideLimited))
    else:
        of.write("\n!!! There are no genes without APN>5 (flag_rna_detected05==0) for any xsome\n")
    if not detectedXsomeWideLimited.empty:
        of.write("\nGenes with APN>5 (flag_rna_detected05==1) by xsome\n{}\n".format(detectedXsomeWideLimited))
    else:
        of.write("\n!!! There are no genes with APN>5 (flag_rna_detected05==1) for any xsome\n")

    ## Sex-bias
    bias = annotDF[(annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")]
    detectedXsomeStackBiased = pd.DataFrame(bias[bias['flag_rna_detected05']==1].groupby(
            ['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
    detectedXsomeWideBiased = detectedXsomeStackBiased.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='xsome', fill_value=0).droplevel(0,axis=1)
    noDetectXsomeStackBiased = pd.DataFrame(bias[bias['flag_rna_detected05']==0].groupby(
            ['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
    noDetectXsomeWideBiased = noDetectXsomeStackBiased.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='xsome', fill_value=0).droplevel(0,axis=1)
    of.write("\n\n#### Sex-biased genes only (gene_rna != male or fem) ####\n\n")
    of.write("\nGenes not APN>5 (flag_rna_detected05==0) by xsome\n{}\n".format(noDetectXsomeWideBiased))
    of.write("\nGenes with APN>5 (flag_rna_detected05==1) by xsome\n{}\n".format(detectedXsomeWideBiased))

    of.close()
  
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()