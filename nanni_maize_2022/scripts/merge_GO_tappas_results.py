#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge tappAS DEA output files with detection flags")

    # Input data
    parser.add_argument("-t", "--tappas", dest="inTappas", required=True, help="Input merged tappAS gene-level result CSV file")
    parser.add_argument("--type", dest="inType", choices=['DE','DIU'], required=False, default="DIU", help="Type of tappAS file (DE or DIU), default is DIU")
    parser.add_argument("--GOid", dest="inID", required=True, help="Input merged GO ID CSV file")
    parser.add_argument("--GOterm", dest="inTerm", required=True, help="Input GO term TSV file")

    # Output data
    parser.add_argument("--output-full", dest="outFull", required=True, help="Output file for full table")
    parser.add_argument("--output-5", dest="outFive", required=False, help="Output file for genes detected and DE in all 5 genotypes")

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

def main():
    # Get input files
    tappasDF = pd.read_csv(args.inTappas,low_memory=False)
    idDF = pd.read_csv(args.inID, low_memory=False)
    termDF = pd.read_csv(args.inTerm, sep="\t", quoting=3, low_memory=False)
    termDF.columns = termDF.columns.str.strip("\"")
    termDF = termDF.apply(lambda x: x.str.strip("\""))
    termDF = termDF[['go_id','Term','Secondary']].drop_duplicates()
    
    # Melt GO dataframe and split GO ID's by "|"
    meltDF = pd.melt(idDF,id_vars="ID",value_vars=['GO_molFunction','GO_bioProcess','GO_cellComponent'])
    splitDF = split_var(df=meltDF,col_name='value',sep="|").rename(columns={
            'ID':'gene_id',
            'value':'GOid',
            'variable':'GOtype'}).dropna()
    
    # Merge with GO terms
    goDF1 = pd.merge(splitDF,termDF,how='outer',left_on='GOid',right_on='go_id',indicator='merge_check1')
#    goDF1['merge_check1'].value_counts()
#    both          635442
#    right_only     36505
#    left_only       5567
    goDF2 = pd.merge(goDF1[goDF1['merge_check1']=='left_only'][splitDF.columns],
                     termDF,how='left',left_on='GOid',right_on='Secondary',indicator='merge_check2')
#    goDF2['merge_check2'].value_counts()
#    both          5524
#    left_only       43
#    right_only       0
#    goDF2[goDF2['merge_check2']=="left_only"]['GOid'].unique()

# 2 ID's not found to have terms ['GO:0003840', 'GO:0000059'], upon searching these are obsolete terms (dropped)
    # Get final go dataframe, make lists of GOid and GOterm separated by "|", then make wide unique on gene_id
    goDF1 = goDF1[goDF1['merge_check1']=='both'].drop(columns=['merge_check1','Secondary','go_id'])
    goDF2 = goDF2[goDF2['merge_check2']=='both'].drop(columns=['merge_check2','Secondary','go_id'])
    goDF = pd.concat([goDF1,goDF2],ignore_index=True,sort=True).rename(columns={'Term':'GOterm'}).drop_duplicates()
    listDF = goDF.groupby(['gene_id','GOtype'])[['GOid','GOterm']].agg(lambda x: "|".join(x)).reset_index()
    listDF['GOtype'] = listDF['GOtype'].str.strip("GO_")
    wideDF = listDF.pivot(index='gene_id',columns='GOtype',values=['GOid','GOterm'])
    wideDF.columns = ["_".join(col).strip() for col in wideDF.columns.values]
    wideDF = wideDF.reset_index()
    
    # Merge results and GO annotations
    mergeDF = pd.merge(wideDF,tappasDF,how='outer',on='gene_id',indicator='merge_check')

    # Output how many genes do not have GO annotation
    print("{} genes with no GO annotation\n{} genes with at least one GO annotation".format(
            mergeDF['merge_check'].value_counts()['right_only'],mergeDF['merge_check'].value_counts()['both']))
    
    # Get only lines with tappas (both and right_only)
    mergeDF = mergeDF[(mergeDF['merge_check']=="both")|(mergeDF['merge_check']=="right_only")].drop(columns=['merge_check'])
    
    # Add flag for all 5 genotype DE
    if args.inType == "DE":
        if ('flag_detect_DE_all5' not in mergeDF.columns) and ('flag_detect_DE_B73' in mergeDF.columns):
            mergeDF['flag_detect_DE_all5'] = np.where((mergeDF['flag_detect_DE_B73']==1)&(mergeDF['flag_detect_DE_C123']==1)&
                    (mergeDF['flag_detect_DE_Hp301']==1)&(mergeDF['flag_detect_DE_Mo17']==1)&
                    (mergeDF['flag_detect_DE_NC338']==1),1,0)
    
    # Output merged file of GO annoations and tappas results
    mergeDF.to_csv(args.outFull,index=False)
    
    # Output merged file of GO annotations and tappas results for
    #     the 83 genes detected and DE in all 5 genotypes
    if args.inType == "DE" and args.outFive is not None:
        mergeDF[mergeDF['flag_detect_DE_all5']==1].to_csv(args.outFive,index=False)
            
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

