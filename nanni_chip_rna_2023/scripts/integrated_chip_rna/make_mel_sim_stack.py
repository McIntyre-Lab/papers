#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get association of gene sex biased expression and the presence of chromatin marks")

    # Input data
    parser.add_argument("-m", "--mel-input", dest="inFileM", required=True, help="CSV file of Dmel gene-level annotation flags with combo flag gene_sex_bias_ttest_foldchange")
    parser.add_argument("-s", "--sim-input", dest="inFileS", required=True, help="CSV file of Dsim gene-level annotation flags with combo flag gene_sex_bias_ttest_foldchange")
    parser.add_argument("--ortholog", dest="inOrtho", required=True, help="CSV file of orthologous gene in Dmel and Dsim")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file for counts")

    args = parser.parse_args()
    return args

def main():
    # Get input files
    annotMel = pd.read_csv(args.inFileM, low_memory=False)
    annotSim = pd.read_csv(args.inFileS, low_memory=False)
    orthoDF = pd.read_csv(args.inOrtho, low_memory=False)

    # Values in gene_sex_bias_ttest_foldchange: male_and_female, male, female, unbiased, ttestLowFC, nan
    # Values in gene_k4_sex: male_and_female, male, fem, unb, none
    # Coordinates are [gene_sex_bias_ttest_foldchange, gene_k4_sex]
    cols = ['fbgn','xsome','gene_rna','gene_k4_sex','gene_k27_sex','gene_trend_ttest']

    # Subset columns and add species variable
    annotMel = annotMel[cols]
    annotSim = annotSim[cols]
    annotMel['species'] = "mel"
    annotSim['species'] = "sim"
    
    # Merge on orthologs to get ortholog flag
    mergeMelOrtho = pd.merge(annotMel,orthoDF['mel_geneID'].drop_duplicates(),how='left',left_on='fbgn',right_on='mel_geneID',indicator='merge_check')
    mergeMelOrtho['flag_ortholog'] = np.where(~mergeMelOrtho['mel_geneID'].isna(),1,0)
    mergeMelOrtho = mergeMelOrtho.drop(columns=['mel_geneID','merge_check'])
    mergeSimOrtho = pd.merge(annotSim,orthoDF['sim_geneID'].drop_duplicates(),how='left',left_on='fbgn',right_on='sim_geneID',indicator='merge_check')
    mergeSimOrtho['flag_ortholog'] = np.where(~mergeSimOrtho['sim_geneID'].isna(),1,0)
    mergeSimOrtho = mergeSimOrtho.drop(columns=['sim_geneID','merge_check'])
    
    # Concat mel and sim
    annotAll = pd.concat([mergeMelOrtho,mergeSimOrtho],ignore_index=True)
    
    # Set flag for having a m/f k4 or k27 mark
    annotAll['flag_has_f_k4'] = np.where(annotAll['gene_k4_sex'].isin(['fem','male_and_female','unb']),1,0)
    annotAll['flag_has_m_k4'] = np.where(annotAll['gene_k4_sex'].isin(['male','male_and_female','unb']),1,0)
    annotAll['flag_has_f_k27'] = np.where(annotAll['gene_k27_sex'].isin(['fem','male_and_female','unb']),1,0)
    annotAll['flag_has_m_k27'] = np.where(annotAll['gene_k27_sex'].isin(['male','male_and_female','unb']),1,0)
    
    # Flag overcompensated and undercompensated genes (using trend_ttest only)
    overConditions = [(annotAll['xsome']=="X")&(annotAll['gene_trend_ttest']=="male")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none")),
                      (annotAll['xsome']=="X")&(annotAll['gene_trend_ttest']!="male")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none"))]
    overChoices = ["1","0"]
    annotAll['flag_overcompensated'] = np.select(overConditions,overChoices,"")
    underConditions = [(annotAll['xsome']=="X")&(annotAll['gene_trend_ttest']=="female")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none")),
                      (annotAll['xsome']=="X")&(annotAll['gene_trend_ttest']!="female")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none"))]
    underChoices = ["1","0"]
    annotAll['flag_undercompensated'] = np.select(underConditions,underChoices,"")
    
    # Get sex-biased genes on the Autosomes
    femConditions = [(annotAll['xsome']=="A")&(annotAll['gene_trend_ttest']=="male")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none")),
                      (annotAll['xsome']=="A")&(annotAll['gene_trend_ttest']!="male")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none"))]
    femChoices = ["1","0"]
    annotAll['flag_female_bias'] = np.select(femConditions,femChoices,"")
    maleConditions = [(annotAll['xsome']=="A")&(annotAll['gene_trend_ttest']=="female")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none")),
                      (annotAll['xsome']=="A")&(annotAll['gene_trend_ttest']!="female")&((annotAll['gene_rna']!="fem")&(annotAll['gene_rna']!="male")&(annotAll['gene_rna']!="none"))]
    maleChoices = ["1","0"]
    annotAll['flag_male_bias'] = np.select(maleConditions,maleChoices,"")

    annotAll.to_csv(args.outFile,index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()