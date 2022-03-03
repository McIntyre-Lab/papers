#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge tappAS DIU output files with detection flags")

    # Input data
    parser.add_argument("-t", "--input-directory", dest="inDir", required=True, help="Input directory of maize tappAS output")
    parser.add_argument("-f", "--flag", dest="inFlag", required=True, help="Input TSV of on/off flags")
    parser.add_argument("-a", "--annot", dest="inAnnot", required=True, help="Input CSV of transcript_id to gene_id pairs in transcriptome")
    parser.add_argument("-e", "--exclude", dest="exclude", required=False, action='append', help="Samples to exclude from expression matrices, multiple values can be listed with each '-e'")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file")
    parser.add_argument("-m", "--major-switch-output", dest="outMajor", required=True, help="Output file of only genes DIU and major isoform switching in at least one genotype")

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
    for genotype in ["C123","Hp301","Mo17","NC338"]:
        tempDF = pd.read_csv("{}/B73_vs_{}_tappAS_DIUGene_Transcripts.tsv".format(args.inDir,genotype),sep="\t")
        tempDF['flag_DIU_'+genotype] = np.where(tempDF['DIU Result']=="DIU",1,0)
        tempDF['flag_DIU_majorIsoformSwitch_'+genotype] = np.where((tempDF['flag_DIU_'+genotype]==1)&
              (tempDF['Major Isoform Switching']=="YES"),1,0)
        tempDF = tempDF.rename(columns={'#Gene':'gene_id',
                                        'Q-Value':'DIU_qval_'+genotype,
                                        'Total Usage Change':'total_usage_change_'+genotype,
                                        'B73 MeanExpLevel':'mean_TPM_ambient_B73',
                                        genotype+' MeanExpLevel':'mean_TPM_ambient_'+genotype})
        tempDF = tempDF[['gene_id','DIU_qval_'+genotype,'flag_DIU_'+genotype,'flag_DIU_majorIsoformSwitch_'+genotype,
                         'total_usage_change_'+genotype,'mean_TPM_ambient_B73','mean_TPM_ambient_'+genotype]]
        tempMerge = pd.merge(mergeDF,tempDF,how='outer',on="gene_id",indicator='merge_check')
#        tempMerge['merge_check'].value_counts()
        print("\t{} genes not detected in {}".format(len(tempMerge[tempMerge['merge_check']=='left_only']),genotype))
        del(tempMerge['merge_check'])
        tempMerge['flag_detect_DIU_majorIsoformSwitch_'+genotype] = np.where((tempMerge['flag_DIU_majorIsoformSwitch_'+genotype]==1)&
                (tempMerge['flag_B73_Amb']+tempMerge['flag_'+genotype+'_Amb']>0),1,0)
        tempMerge['flag_both_detect_DIU_majorIsoformSwitch_'+genotype] = np.where((tempMerge['flag_DIU_majorIsoformSwitch_'+genotype]==1)&
                (tempMerge['flag_B73_Amb']+tempMerge['flag_'+genotype+'_Amb']==2),1,0)
        mergeDF = tempMerge.copy()
    
    # Output merged file of flags and tappas results
    mergeDF.to_csv(args.outFile,index=False)
    mergeDF['sum_detect_DIU_majorIsoformSwitch'] = mergeDF[[c for c in mergeDF.columns if "flag_detect_DIU_majorIsoformSwitch" in c]].sum(axis=1)
    mergeDF['sum_both_detect_DIU_majorIsoformSwitch'] = mergeDF[[c for c in mergeDF.columns if "flag_both_detect_DIU_majorIsoformSwitch" in c]].sum(axis=1)
    
    # Print out counts
    print("\nNumber of DIU genes with major isoform switch detected between ambient B73 and ambient of each genotype:\n{}".format(
            mergeDF[[c for c in mergeDF.columns if "flag_detect_DIU_majorIsoformSwitch" in c]].sum().to_string()))
    print("\nNumber of genes shared among genotypes:\nnum_geno   num_gene\n{}".format(
            mergeDF['sum_detect_DIU_majorIsoformSwitch'].value_counts().to_string()))
    # The 1 gene with 3 is Zm00001d054001 with Hp301, Mo17, and NC338
    
    # Output file of only genes DIU and major isoform switching in at least one genotype
    mergeDF[mergeDF['sum_detect_DIU_majorIsoformSwitch']>0].to_csv(args.outMajor,index=False)
            
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

