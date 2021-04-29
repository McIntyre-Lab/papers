#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Plot boxplot of gene mean expression TPM values (sum of transcript means) within the groups of genes 1) shared by Wang 2018 leaf transcriptome and 5 genotype/2 treatment leaf transcriptome, and 2) only in 5 genotype/2 treatment leaf transcriptome")

    # Input data
    parser.add_argument("-f", "--flag", dest="inFlag", required=True, help="Input TSV of combined rsem expression matrix with TPM on/off flag")
    parser.add_argument("-a", "--annot", dest="inAnnot", required=True, help="Input CSV of transcript_id to gene_id pairs in transcriptome")
    parser.add_argument("-1", "--input-file1", dest="inF1", required=True, help="CSV file for 5 genotype/2 treatment leaf transcrioptome of event_id to transcript_id to gene_id (*_event2transcript2gene_index.csv)")
    parser.add_argument("-2", "--input-file2", dest="inF2", required=True, help="CSV file for Wang 2018 B73 leaf transcriptome of event_id to transcript_id to gene_id (*_event2transcript2gene_index.csv)")
    parser.add_argument("--missing1", dest="inM1", required=True, help="List of genes missing from EA annotations in 5 genotype/2 treatment leaf transcrioptome")
    parser.add_argument("--missing2", dest="inM2", required=True, help="List of genes missing from EA annotations in Wang 2018 B73 leaf transcriptome")

    # Output data
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Output prefix")

    args = parser.parse_args()
    return args

def main():
    # Get on/off flags, transcript_id to gene_id pairs,
    #     and event_id to transcript_id to gene_id files from 5 genotype (compareDF1)
    #     and Wang 2018 leaf dataset (compareDF2)
    flagDF = pd.read_csv(args.inFlag, sep="\t",low_memory=False)
    annotDF = pd.read_csv(args.inAnnot)
    compareDF1 = pd.read_csv(args.inF1, low_memory=False)
    compareDF2 = pd.read_csv(args.inF2, low_memory=False)
    missingDF1 = pd.read_csv(args.inM1,names=['gene_id'])
    missingDF1['transcript_id'] = 1
    missingDF2 = pd.read_csv(args.inM2,names=['gene_id'])
    missingDF2['transcript_id'] = 1

    # Extract all unique pairs of individual transcript_id to gene_id
    #     and count transcripts per gene in each dataset
    xcrptDF1 = compareDF1[(~compareDF1['transcript_id'].str.contains("|",regex=False))&
            (~compareDF1['gene_id'].str.contains("|",regex=False))&
            (compareDF1['transcript_id']!="Unannotated")][['gene_id','transcript_id']].drop_duplicates()
    geneDF1 = xcrptDF1.groupby('gene_id')['transcript_id'].count().reset_index()
    geneDF1 = pd.concat([geneDF1,missingDF1])
    xcrptDF2 = compareDF2[(~compareDF2['transcript_id'].str.contains("|",regex=False))&
            (~compareDF2['gene_id'].str.contains("|",regex=False))&
            (compareDF2['transcript_id']!="Unannotated")][['gene_id','transcript_id']].drop_duplicates()
    geneDF2 = xcrptDF2.groupby('gene_id')['transcript_id'].count().reset_index()
    geneDF2 = pd.concat([geneDF2,missingDF2])
    
    # Merge datasets and count
    mergeDF = pd.merge(geneDF1,geneDF2,how='outer',on='gene_id',suffixes=('_1','_2'),indicator='merge_check')
    mergeDF['transcript_id_2'] = np.where(mergeDF['merge_check']=='left_only',0,mergeDF['transcript_id_2'])
    countDF = mergeDF['merge_check'].value_counts()
    print("{} genes in 5 genotype/2 treatment leaf transcrioptome only\n{} genes in Wang 2018 B73 leaf transcriptome only\n{} genes in both".format(
            countDF['left_only'],countDF['right_only'],countDF['both']))
    
    # Merge gene_id into flags
    geneDF = pd.merge(annotDF,flagDF,how='outer',on='transcript_id',indicator='merge_check')
#    geneDF['merge_check'].value_counts()
    # Fix the 89 single transcript genes missing in EA files
    #     (have the same value for transcript_id and gene_id)
    geneDF = geneDF[geneDF['merge_check']!='left_only']
    geneDF['gene_id'] = np.where(geneDF['merge_check']=='right_only',geneDF['transcript_id'],geneDF['gene_id'])
    del(geneDF['merge_check'])
#    geneDF['gene_id'].isna().any()
    
    # Get genes that are detected in both treatments of at least one genotype
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        geneDF['flag_detect_'+genotype] = np.where(geneDF['flag_'+genotype+'_Amb']+
                geneDF['flag_'+genotype+'_Ele']>0,1,0)
        
    # Get max of flags within genes
    geneCollapseDF = geneDF.groupby('gene_id').agg(dict([(c,'max') for c in geneDF.columns if "flag_" in c]+[(c,'sum') for c in geneDF.columns if "mean" in c])).reset_index()

    # Merge flags with dataset counts both or in 5 genotype/2 treatment only
    subsetDF = mergeDF[mergeDF['merge_check']!="right_only"].copy()
    subsetDF['group'] = np.where(subsetDF['merge_check']=="both","Both Transcriptomes","5 Genotypes\nTranscriptome Only")
    del(subsetDF['merge_check'])
    allDF = pd.merge(subsetDF,geneCollapseDF,how='outer',on='gene_id',indicator='merge_check')

    # Mean expression distribution
    meanDF = allDF[[c for c in allDF.columns if ("flag_detect" in c) or ("mean" in c) or (c=="group")]].copy()
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        # Set detected mean to NA if not detected
        meanDF['detect_mean_'+genotype+"_Amb"] = np.where(meanDF['flag_detect_'+genotype]==1,meanDF['mean_'+genotype+'_Amb'],np.nan)
        meanDF['detect_mean_'+genotype+"_Ele"] = np.where(meanDF['flag_detect_'+genotype]==1,meanDF['mean_'+genotype+'_Ele'],np.nan)
    stackDF = meanDF[[c for c in meanDF.columns if ("detect_mean" in c) or ('group' in c)]].set_index('group').stack().reset_index().rename(columns={'level_1':'sample',0:'detect_mean_TPM'})
    stackDF['sample'] = stackDF['sample'].str.split("_").str[-2:].str.join(" ")
    fig = plt.figure(figsize=(15,7))
    fig.text(0.08, 0.5, "Detected Genotype Mean Gene Expression (TPM)", ha='center', va='center', rotation='vertical')
#    axUpper = plt.subplot2grid((2,1),(0,0),fig=fig)
#    axLower = plt.subplot2grid((2,1),(1,0),fig=fig)
    ax = plt.subplot2grid((1,1),(0,0),fig=fig)
#    axUpper.set_ylim(66,375)
#    axLower.set_ylim(0,65)
    ax.set_ylim(0,375)
    # Color codes from custom palette generated here: https://davidmathlogic.com/colorblind/#%23949494-%23D8D8D8-%23FF8700-%23FFDBB0-%232EBEDC-%23A8F4FF-%2397EC00-%23BDFF7D-%23D836C2-%23D898CD-%23FF0043
    colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
    sns.boxplot(x='group',y='detect_mean_TPM',data=stackDF,hue='sample',fliersize=0,ax=ax,palette=colors)
#    sns.boxplot(x='group',y='detect_mean_TPM',data=stackDF,hue='sample',fliersize=0,ax=axUpper,palette=colors)
#    sns.boxplot(x='group',y='detect_mean_TPM',data=stackDF,hue='sample',fliersize=0,ax=axLower,palette=colors)
#    axUpper.set_xlabel("")
#    axUpper.set_xticks([])
#    axUpper.set_ylabel("")
#    axUpper.legend(title="")
    ax.legend(title="")
#    axLower.set_xlabel("Number of Detected Genotypes")
    ax.set_xlabel("")
#    axLower.set_ylabel("")
    ax.set_ylabel("")
#    axLower.legend().set_visible(False)
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("{}_gene_meanExp_boxplot.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_gene_meanExp_boxplot.pdf".format(args.outPrefix),dpi=600,format="pdf")

    # Merge in names of genes to check lists
#    geneNameDF = pd.merge(nameDF,geneCollapseDF,how='outer',on='gene_id',indicator='merge_check')
#    geneNameDF = geneNameDF[geneNameDF['merge_check']=='both']
#    del(geneNameDF['merge_check'])
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
