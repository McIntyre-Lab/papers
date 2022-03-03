#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Plot density for expression ratio of IR transcripts vs. total expression within the gene in each genotype")

    # Input data
    parser.add_argument("-e", "--expression", dest="inExp", required=True, help="Input TSV of combined rsem isoform expression matrix with TPM on/off flag")
    parser.add_argument("-d", "--distance", dest="inDist", required=True, help="Input CSV of pairwise transcript distances (*_pairwise_transcript_distance.csv)")
    parser.add_argument("-f", "--fragment", dest="inFrag", required=True, help="Input CSV of intron retention flagged fragment annotation file (*_exon_fragment_annotations_flag_IR.csv)")

    # Output data
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Output prefix")

    args = parser.parse_args()
    return args

def split_var(df,col_name=None,sep=None,sort_list=None):
    # Split variable by sep and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def unnest(df, explode):
    idx = df.index.repeat(df[explode[0]].str.len())
    df1 = pd.concat([
        pd.DataFrame({x: np.concatenate(df[x].values)}) for x in explode], axis=1)
    df1.index = idx
    return df1.join(df.drop(explode, 1), how='left')

def main():
    # Get input pairwise distance file and check if it is empty
    distDF = pd.read_csv(args.inDist, low_memory=False)
    if len(distDF) == 0:
        print("ERROR: Empty distance file provided ... No plots generated")
        sys.exit()
        
    # Get number of transcripts per gene from distance file
    # Number of pairwise comparisons p is p = k(k-1)/2 where k is the number of transcripts
    #     so k = (1 + sqrt(1+8p))/2
    distDF['numPair'] = distDF.groupby('gene_id')['transcript_1'].transform('count')
    distDF['numXcrpt'] = (1+np.sqrt(1+(8*distDF['numPair'])))/2
    del(distDF['numPair'])

    # Get input fragment annotation file flagged with involvement in IR
    fragDF = pd.read_csv(args.inFrag, low_memory=False)
    
    # Flag NIC/NNC transcripts ("PB" in transcript_id)
    fragDF['flag_has_NIC_NNC'] = np.where(fragDF['transcript_id'].str.contains("PB"),1,0)
    
    # Flag if NIC/NNC has putative IR fusion
    fragDF['flag_NIC_NNC_IR_fusion'] = np.where((fragDF['flag_has_NIC_NNC']==1)&
          (fragDF['flag_intron_retention']==1),1,0)
    
    # Split fragments by transcript_id
    splitDF = split_var(df=fragDF,col_name="transcript_id",sep="|")
    # Fix IR flag to only flag the transcript with the IR event (will really only work for 2 transcript genes with no multigene)
    splitDF['flag_IR_isoform'] = np.where((splitDF['flag_intron_retention']==1)&
                       (splitDF['annotation_frequency']=="Unique"),1,0)
    
    # Get unique transcript_id values with IR flags
    xcrptFragDF = splitDF[splitDF['flag_multigene']==0].groupby('transcript_id')[['flag_intron_retention',
                               'flag_has_NIC_NNC','flag_NIC_NNC_IR_fusion',
                               'flag_IR_isoform']].max().reset_index()

    # Flag alternative exons (prop_ER_diff > 0)
    distDF['flag_alt_exon'] = np.where(distDF['prop_ER_diff']>0,1,0)
    
    # Flag alternative donor/acceptors (num_nt_diff_in_shared_ER > 0)
    distDF['flag_alt_donor_acceptor'] = np.where(distDF['num_nt_diff_in_shared_ER']>0,1,0)
     
    # Get unique transcript_id values with distance flags
    distDF['transcript_id'] = distDF[['transcript_1','transcript_2']].values.tolist()
    xcrptDistDF = unnest(distDF,['transcript_id']).groupby('transcript_id')[['gene_id','flag_alt_exon','flag_alt_donor_acceptor',
                               'numXcrpt','num_nt_diff_in_shared_ER']].max().reset_index().rename(columns={'num_nt_diff_in_shared_ER':'max_num_nt_diff_in_shared_ER'})
    
    # Merge transcript dataframes on transcript_id and only keep those in distance file (have >=2 transcripts)
    mergeDF = pd.merge(xcrptFragDF,xcrptDistDF,how='outer',on='transcript_id',indicator='merge_check')
    mergeDF = mergeDF[mergeDF['merge_check']=="both"]
    del(mergeDF['merge_check'])
    del(xcrptDistDF)
    del(xcrptFragDF)
    
    # Check that no genes have none (all should have at least one type)
    if len(mergeDF[(mergeDF['flag_intron_retention']==0)&(mergeDF['flag_alt_exon']==0)&(mergeDF['flag_alt_donor_acceptor']==0)]) > 0:
        print("ERROR: At least one gene with no splicing classifications found")
        sys.exit()
    
    # Output csv file of flags
#    mergeDF.to_csv(args.outFile,index=False)
    
    
#### Expression ####
    
    # Get input file of isoforms expression and flags
    expDF = pd.read_csv(args.inExp, sep="\t",low_memory=False)
    
    # Merge with splicing flags
    expFlagDF = pd.merge(mergeDF,expDF,how='outer',on='transcript_id',indicator='merge_check')

    # Get ratio of gene expression for each transcript within a sample
    # Get transcripts that are detected in both treatments of at least one genotype
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        expFlagDF['total_gene_mean_'+genotype+'_Amb'] = expFlagDF.groupby('gene_id')['mean_'+genotype+'_Amb'].transform('sum')
        expFlagDF['ratio_mean_'+genotype+'_Amb'] = expFlagDF['mean_'+genotype+'_Amb']/expFlagDF['total_gene_mean_'+genotype+'_Amb']
        expFlagDF['total_gene_mean_'+genotype+'_Ele'] = expFlagDF.groupby('gene_id')['mean_'+genotype+'_Ele'].transform('sum')
        expFlagDF['ratio_mean_'+genotype+'_Ele'] = expFlagDF['mean_'+genotype+'_Ele']/expFlagDF['total_gene_mean_'+genotype+'_Ele']
        expFlagDF['flag_detect_'+genotype] = np.where(expFlagDF['flag_'+genotype+'_Amb']+
                expFlagDF['flag_'+genotype+'_Ele']>0,1,0)
        # Set detected mean to NA if not detected
        expFlagDF['detect_mean_'+genotype+"_Amb"] = np.where(expFlagDF['flag_detect_'+genotype]==1,expFlagDF['mean_'+genotype+'_Amb'],np.nan)
        expFlagDF['detect_ratio_mean_'+genotype+"_Amb"] = np.where(expFlagDF['flag_detect_'+genotype]==1,expFlagDF['ratio_mean_'+genotype+'_Amb'],np.nan)
        expFlagDF['detect_mean_'+genotype+"_Ele"] = np.where(expFlagDF['flag_detect_'+genotype]==1,expFlagDF['mean_'+genotype+'_Ele'],np.nan)
        expFlagDF['detect_ratio_mean_'+genotype+"_Ele"] = np.where(expFlagDF['flag_detect_'+genotype]==1,expFlagDF['ratio_mean_'+genotype+'_Ele'],np.nan)

    # Select transcripts with IR in genes with 2 transcripts (and no alt exon)
    twoXcrptIR = expFlagDF[(expFlagDF['numXcrpt']==2)&(expFlagDF['flag_alt_exon']==0)&
                           (expFlagDF['flag_IR_isoform']==1)].copy()
    print("{} ({} NIC/NNC) genes out of {} ({} NIC/NNC) with two transcripts and IR event (and not alt. exon) have different IR in both transcripts".format(
        twoXcrptIR['gene_id'].duplicated().sum(),
        twoXcrptIR[twoXcrptIR['transcript_id'].str.contains("PB")]['gene_id'].duplicated().sum(),
        twoXcrptIR['gene_id'].nunique(),
        twoXcrptIR[twoXcrptIR['transcript_id'].str.contains("PB")]['gene_id'].nunique()))

    # Mean expression distribution - ALL IR (regardless of detection)
    meanDF = twoXcrptIR[[c for c in twoXcrptIR.columns if ("mean_" in c) and ("detect" not in c) and ("total" not in c) and ("ratio" not in c)]]
    meanDF.columns = meanDF.columns.str.split("_").str[-2:].str.join(" ")
    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot2grid((1,1),(0,0),fig=fig)
    ax.set_xlim(0,50)
    ax.set_xlabel("IR Isoform Mean Expression")
    colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
    meanDF.plot.kde(bw_method=0.01,ax=ax,color=colors)
    plt.title("Mean expression distribution")
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("{}_2xcrptGene_all_IR_isoform_meanExp_density.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_2xcrptGene_all_IR_isoform_meanExp_density.pdf".format(args.outPrefix),dpi=600,format="pdf")
    print("Mean expression distribution - ALL IR (regardless of detection)\n{}".format(
        meanDF.describe().to_string()))

    # Mean expression distribution - ALL IR (detected only)
    meanDetectDF = twoXcrptIR[[c for c in twoXcrptIR.columns if "detect_mean_" in c]]
    print("{} transcripts detected in at least one treatment out of {} transcripts with IR in genes with 2 transcripts total".format(
        (~meanDetectDF.isna()).all(axis=1).sum(),len(meanDetectDF)))
    meanDetectDF.columns = meanDetectDF.columns.str.split("_").str[-2:].str.join(" ")
    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot2grid((1,1),(0,0),fig=fig)
    ax.set_xlim(0,50)
    ax.set_xlabel("IR Isoform Mean Expression")
    colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
    meanDetectDF.plot.kde(bw_method=0.01,ax=ax,color=colors)
    plt.title("Mean expression distribution")
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("{}_2xcrptGene_detected_IR_isoform_meanExp_density.png".format(args.outPrefix),dpi=600,format="png")
#    plt.savefig("{}_detected_IR_isoform_meanExp_density.pdf".format(args.outPrefix),dpi=600,format="pdf")
    print("Mean expression distribution - ALL IR (detected only)\n{}".format(
        meanDetectDF.describe().to_string()))
    
    # Ratio of IR expression distribution - ALL IR (regardless of detection)
    ratioDF = twoXcrptIR[[c for c in twoXcrptIR.columns if ("ratio_" in c) and ("detect" not in c)]]
    ratioDF.columns = ratioDF.columns.str.split("_").str[-2:].str.join(" ")
    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot2grid((1,1),(0,0),fig=fig)
    ax.set_xlim(0,1)
    ax.set_xlabel("IR/total Expression Ratio")
    colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
    ratioDF.plot.kde(ax=ax,color=colors)
    plt.title("Ratio of expression distribution")
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("{}_2xcrptGene_all_IR_isoform_ratioExp_density.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_2xcrptGene_all_IR_isoform_ratioExp_density.pdf".format(args.outPrefix),dpi=600,format="pdf")
    print("Ratio of expression distribution - ALL IR (regardless of detection)\n{}".format(
        ratioDF.describe().to_string()))

    # Ratio of expression distribution - ALL IR (detected only)
    ratioDetectDF = twoXcrptIR[[c for c in twoXcrptIR.columns if "detect_ratio_" in c]]
    ratioDetectDF.columns = ratioDetectDF.columns.str.split("_").str[-2:].str.join(" ")
    fig = plt.figure(figsize=(15,7))
    ax = plt.subplot2grid((1,1),(0,0),fig=fig)
    ax.set_xlim(0,1)
    ax.set_xlabel("IR/total Expression Ratio")
    colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
    ratioDetectDF.plot.kde(ax=ax,color=colors)
    plt.title("Ratio of expression distribution")
    plt.subplots_adjust(hspace=0.05)
    plt.savefig("{}_2xcrptGene_detected_IR_isoform_ratioExp_density.png".format(args.outPrefix),dpi=600,format="png")
#    plt.savefig("{}_detected_IR_isoform_meanExp_density.pdf".format(args.outPrefix),dpi=600,format="pdf")
    print("Ratio of expression distribution - ALL IR (detected only)\n{}".format(
        ratioDetectDF.describe().to_string()))

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
