#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from upsetplot import UpSet
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Plot presence/absence of genes in maize PacBio transcriptome")

    # Input data
    parser.add_argument("-f", "--flag", dest="inFlag", required=True, help="Input TSV of combined rsem expression matrix with TPM on/off flag")
    parser.add_argument("-a", "--annot", dest="inAnnot", required=True, help="Input CSV of transcript_id to gene_id pairs in transcriptome")
    parser.add_argument("-t", "--tappas", dest="inTappas", required=True, help="Input CSV of tappas gene-level DE results")
    parser.add_argument("-n", "--name", dest="inName", required=True, help="Input TSV of gene_id (column 1) and gene_name (column 2) without header")
    parser.add_argument("-v", "--p-val", dest="inPval", action='store_true', required=False, help="Output histogram plots of DE result FDR corrected p-values for each genotype")
    parser.add_argument("-m", "--mean-exp", dest="inMeanExp", action='store_true', required=False, help="Output boxplot of gene mean expression TPM values (sum of transcript means) within the groups of number of genes detected (only plot expression values of detected genotype - aka at least one treatment detected)")

    # Output data
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Output prefix")

    args = parser.parse_args()
    return args

def main():
    # Get on/off flags, tappas results, transcript_id to gene_id pairs,
    #     and gene_id to gene_name pairs
    flagDF = pd.read_csv(args.inFlag, sep="\t",low_memory=False)
    tappasDF = pd.read_csv(args.inTappas,low_memory=False)
    annotDF = pd.read_csv(args.inAnnot)
    nameDF = pd.read_csv(args.inName,sep="\t",names=['gene_id','gene_name'])
    
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
    # Make extra flags and count
    # Amb vs Ele
#    geneCollapseDF = geneDF.groupby('gene_id')[[c for c in geneDF.columns if "flag_" in c]].max().reset_index()
    geneCollapseDF = geneDF.groupby('gene_id').agg(dict([(c,'max') for c in geneDF.columns if "flag_" in c]+[(c,'sum') for c in geneDF.columns if "mean" in c])).reset_index()
    geneCollapseDF['flag_all_on_Amb'] = np.where(geneCollapseDF[[c for c in geneCollapseDF.columns if ('flag_' in c) and ('_Amb'in c)]].sum(axis=1)==5,1,0)
    geneCollapseDF['flag_all_on_Ele'] = np.where(geneCollapseDF[[c for c in geneCollapseDF.columns if ('flag_' in c) and ('_Ele'in c)]].sum(axis=1)==5,1,0)
    geneCollapseDF['flag_all_off_Amb'] = np.where(geneCollapseDF[[c for c in geneCollapseDF.columns if ('flag_' in c) and ('_Amb'in c) and ("all_on" not in c)]].sum(axis=1)==0,1,0)
    geneCollapseDF['flag_all_off_Ele'] = np.where(geneCollapseDF[[c for c in geneCollapseDF.columns if ('flag_' in c) and ('_Ele'in c) and ("all_on" not in c)]].sum(axis=1)==0,1,0)
    print(pd.crosstab(geneCollapseDF['flag_all_on_Amb'],geneCollapseDF['flag_all_off_Ele']))
    print("\n{}".format(pd.crosstab(geneCollapseDF['flag_all_off_Amb'],geneCollapseDF['flag_all_on_Ele'])))

    # Merge in names of genes to check lists
    geneNameDF = pd.merge(nameDF,geneCollapseDF,how='outer',on='gene_id',indicator='merge_check')
    geneNameDF = geneNameDF[geneNameDF['merge_check']=='both']
    del(geneNameDF['merge_check'])

    # Output expression matrices for genes detected in only Amb and only Ele
    geneNameDF[(geneNameDF['flag_all_on_Amb']==1)&
                   (geneNameDF['flag_all_off_Ele']==1)].to_csv("{}_Amb_only.csv".format(args.outPrefix),index=False)
    geneNameDF[(geneNameDF['flag_all_on_Ele']==1)&
                   (geneNameDF['flag_all_off_Amb']==1)].to_csv("{}_Ele_only.csv".format(args.outPrefix),index=False)
    
    # Genotypes
    geneCollapseDF['sum_detect'] = geneCollapseDF[[c for c in geneCollapseDF.columns if "flag_detect" in c]].sum(axis=1)
    print("\nnum_geno  num_gene\n{}".format(geneCollapseDF['sum_detect'].value_counts(sort=False).to_string()))
    print("\nCount genes detected in single genotype:\n{}".format(geneCollapseDF[geneCollapseDF['sum_detect']==1][[c for c in geneCollapseDF.columns if "flag_detect" in c]].sum().to_string()))

    # Mean expression distribution
    if args.inMeanExp:
        meanDF = geneCollapseDF[[c for c in geneCollapseDF.columns if ("flag_detect" in c) or ("mean" in c) or (c=="sum_detect")]].copy()
        for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
            # Set detected mean to NA if not detected
            meanDF['detect_mean_'+genotype+"_Amb"] = np.where(meanDF['flag_detect_'+genotype]==1,meanDF['mean_'+genotype+'_Amb'],np.nan)
            meanDF['detect_mean_'+genotype+"_Ele"] = np.where(meanDF['flag_detect_'+genotype]==1,meanDF['mean_'+genotype+'_Ele'],np.nan)
        stackDF = meanDF[[c for c in meanDF.columns if ("detect_mean" in c) or ('sum_detect' in c)]].set_index('sum_detect').stack().reset_index().rename(columns={'level_1':'sample',0:'detect_mean_TPM'})
        stackDF['sample'] = stackDF['sample'].str.split("_").str[-2:].str.join(" ")
        fig = plt.figure(figsize=(15,7))
        fig.text(0.08, 0.5, "Distribution of Mean Gene Expression (TPM)", ha='center', va='center', rotation='vertical')
        axUpper = plt.subplot2grid((2,1),(0,0),fig=fig)
        axLower = plt.subplot2grid((2,1),(1,0),fig=fig)
        axUpper.set_ylim(66,375)
        axLower.set_ylim(0,65)
        # Color codes from custom palette generated here: https://davidmathlogic.com/colorblind/#%23949494-%23D8D8D8-%23FF8700-%23FFDBB0-%232EBEDC-%23A8F4FF-%2397EC00-%23BDFF7D-%23D836C2-%23D898CD-%23FF0043
        colors = ["#2EBEDC","#A8F4FF","#D836C2","#D898CD","#949494","#D8D8D8","#97EC00","#BDFF7D","#FF8700","#FFDBB0"]
        sns.boxplot(x='sum_detect',y='detect_mean_TPM',data=stackDF,hue='sample',fliersize=0,ax=axUpper,palette=colors)
        sns.boxplot(x='sum_detect',y='detect_mean_TPM',data=stackDF,hue='sample',fliersize=0,ax=axLower,palette=colors)
        axUpper.set_xlabel("")
        axUpper.set_xticks([])
        axUpper.set_ylabel("")
        axUpper.legend(title="")
        axLower.set_xlabel("Number of Genotypes")
        axLower.set_ylabel("")
        axLower.legend().set_visible(False)
        axLower.set_xticklabels(['1\nn='+str(geneCollapseDF['sum_detect'].value_counts()[1]),
                                 '2\nn='+str(geneCollapseDF['sum_detect'].value_counts()[2]),
                                 '3\nn='+str(geneCollapseDF['sum_detect'].value_counts()[3]),
                                 '4\nn='+str(geneCollapseDF['sum_detect'].value_counts()[4]),
                                 '5\nn='+str(geneCollapseDF['sum_detect'].value_counts()[5])])
        plt.subplots_adjust(hspace=0.05)
        plt.savefig("{}_genotype_detect_gene_meanExp_boxplot.png".format(args.outPrefix),dpi=600,format="png")
        plt.savefig("{}_genotype_detect_gene_meanExp_boxplot.pdf".format(args.outPrefix),dpi=600,format="pdf")
    
    # DE plot
    # Requiring that genes detected in at least one treatment of the genotype 
    #     to be considered DE
    detectDF = tappasDF.copy()
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        detectDF['flag_detect_DE_'+genotype] = np.where((detectDF['flag_DE_'+genotype]==1)&
                (detectDF['flag_'+genotype+'_Amb']+detectDF['flag_'+genotype+'_Ele']>0),True,False)
    detectCountDF = detectDF[[c for c in detectDF.columns if "flag_detect_DE" in c]].copy()
    detectCountDF.columns = detectCountDF.columns.str.split("_").str[3]
    detectConcatDF = pd.concat([detectDF,detectCountDF],axis=1)
    detectConcatDF = detectConcatDF.set_index(list(detectCountDF.columns))
    upset = UpSet(detectConcatDF,subset_size='count',show_counts=True,sort_by='degree')
    upset.plot()
    plt.savefig("{}_DE_atLeast1_trt_detected_upset_plot.png".format(args.outPrefix),dpi=600,format="png")
    upset = UpSet(detectConcatDF,subset_size='count',show_counts=True,sort_by='cardinality')
    upset.plot()
    plt.savefig("{}_DE_atLeast1_trt_detected_upset_plot_sortVal.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_DE_atLeast1_trt_detected_upset_plot_sortVal.pdf".format(args.outPrefix),dpi=600,format="pdf")
    
    # Output histogram plots of DE result FDR corrected p-values for each genotype
    if args.inPval:
        for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
            fig = plt.figure(figsize=(12,10))
            allPval = plt.subplot2grid((5,1),(0,0),rowspan=2,fig=fig)
            detectPval = plt.subplot2grid((5,1),(3,0),rowspan=2,fig=fig)
            allPval.hist(detectDF[~detectDF['DE_pval_'+genotype].isna()]['DE_pval_'+genotype],bins=200)
            allPval.set_title("{} DE p-value (FDR) for all genes tested".format(genotype))
            allPval.set_xlabel("p-value")
            allPval.set_ylabel("# of genes")
            detectPval.hist(detectDF[detectDF['flag_detect_DE_'+genotype]==1]['DE_pval_'+genotype],bins=200)
            detectPval.set_title("{} DE p-value (FDR) for all genes detected\nin at least one treatment".format(genotype))
            detectPval.set_xlabel("p-value")
            detectPval.set_ylabel("# of genes")   
            plt.savefig("{}_{}_DE_pval_hist.png".format(args.outPrefix,genotype),dpi=600,format="png")
                           
    # Merge list of DE in all 5 genotypes with gene names
    listMerge = pd.merge(nameDF,detectConcatDF.loc[(True,True,True,True,True),:],
                                                     how='outer',on='gene_id',indicator='merge_check')
    #listMerge['merge_check'].value_counts()
    # merge is good
    
    # Output list
    listMerge[listMerge['merge_check']=="both"][['gene_id','gene_name']+
              [c for c in listMerge.columns if ("Log2FC" in c)or("mean_TPM" in c)]].to_csv(
              "{}_DE_all5genotype.csv".format(args.outPrefix),index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
