#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge tappAS DEA output files with detection flags")

    # Input data
    parser.add_argument("-t", "--input-directory", dest="inDir", required=True, help="Input directory of maize tappAS output")
    parser.add_argument("-f", "--flag", dest="inFlag", required=True, help="Input TSV of on/off flags")
    parser.add_argument("-a", "--annot", dest="inAnnot", required=True, help="Input CSV of transcript_id to gene_id pairs in transcriptome")
    parser.add_argument("-e", "--exclude", dest="exclude", required=False, action='append', help="Samples to exclude from expression matrices, multiple values can be listed with each '-e'")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file")

    args = parser.parse_args()
    return args

def main():
    # Get on/off flags and transcript_id to gene_id pairs
    flagDF = pd.read_csv(args.inFlag, sep="\t",low_memory=False)
    annotDF = pd.read_csv(args.inAnnot)
    
    # Merge gene_id into flags
    geneDF = pd.merge(annotDF,flagDF,how='outer',on='transcript_id',indicator='merge_check')
#    geneDF['merge_check'].value_counts()
    # Fix the 89 single transcript genes missing in EA files
    #     (have the same value for transcript_id and gene_id)
    geneDF = geneDF[geneDF['merge_check']!='left_only']
    geneDF['gene_id'] = np.where(geneDF['merge_check']=='right_only',geneDF['transcript_id'],geneDF['gene_id'])
    del(geneDF['merge_check'])
#    geneDF['gene_id'].isna().any()
    geneDF = geneDF.groupby('gene_id')[[c for c in geneDF.columns if 'flag' in c]].max().reset_index()
    
    # Remove excluded columns if provided
    if args.exclude is not None:
        for s in args.exclude:
            geneDF = geneDF.drop(columns=[c for c in geneDF.columns if c==s])   
            print("Removed {} columns from matrix...".format(s))

    # Count and drop transcripts that are not detected in any samples
    #     and transcripts only detected in one samples
    if 'sum_flag' not in geneDF.columns:
        geneDF['sum_flag'] = geneDF[[c for c in geneDF.columns if "flag_" in c]].sum(axis=1).astype(int)
    print("Detection of genes by number of genotypes:\nnumSamples\tfreq\n{}".format(
            geneDF['sum_flag'].value_counts().sort_index()))
    print("{} transcripts detected in 0 samples or 1 sample only\n...Dropped from tappas files".format(
            len(geneDF[geneDF['sum_flag']<=1])))
    detectDF = geneDF[geneDF['sum_flag']>1].copy().reset_index(drop=True)
    
    # Merge genotype tappAS output files        
    mergeDF = detectDF.copy()
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        tempDF = pd.read_csv("{}/{}_tappAS_DEA_Genes.tsv".format(args.inDir,genotype),sep="\t")
        tempDF['flag_DE_'+genotype] = np.where(tempDF['DEA Result']=="DE",1,0)
        tempDF = tempDF.rename(columns={'#Gene':'gene_id',
                                        '(1 - Probability)':'DE_pval_'+genotype,
                                        'Log2FC':'Log2FC_'+genotype,
                                        'Ambient MeanExpLevel':'mean_TPM_ambient_'+genotype,
                                        'Ozone MeanExpLevel':'mean_TPM_ozone_'+genotype})
        tempDF = tempDF[['gene_id','DE_pval_'+genotype,'flag_DE_'+genotype,
                         'Log2FC_'+genotype,'mean_TPM_ambient_'+genotype,'mean_TPM_ozone_'+genotype]]
        tempMerge = pd.merge(mergeDF,tempDF,how='outer',on="gene_id",indicator='merge_check')
#        tempMerge['merge_check'].value_counts()
        print("\t{} genes not detected in {}".format(len(tempMerge[tempMerge['merge_check']=='left_only']),genotype))
        del(tempMerge['merge_check'])
        tempMerge['flag_detect_DE_'+genotype] = np.where((tempMerge['flag_DE_'+genotype]==1)&
                (tempMerge['flag_'+genotype+'_Amb']+tempMerge['flag_'+genotype+'_Ele']>0),1,0)
        mergeDF = tempMerge.copy()
    
    # Output merged file of flags and tappas results
    mergeDF.to_csv(args.outFile,index=False)
            
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

