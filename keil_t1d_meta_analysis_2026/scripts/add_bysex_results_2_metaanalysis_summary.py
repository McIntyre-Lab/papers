#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 16:25:10 2025

@author: nkeil
"""

import os 
import pandas as pd

ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts"

CD4_DF=pd.read_csv(os.path.join(ind, "metaanalysis_summary_CD4_w_ES_diff.csv"))
CD8_DF=pd.read_csv(os.path.join(ind, "metaanalysis_summary_CD8_w_ES_diff.csv"))

sex_pval_DF=pd.read_csv(os.path.join(ind, "meta_analysis_DE_pdiffs_summary_bysex/meta_analysis_pvals_heterogenity_bysex.csv"))

#Pivot table
sex_pval_DF= sex_pval_DF.pivot_table(
    index=['geneID', 'cellType'],
    columns='sex',
    values=['pval_heterogeneity', 'flag_sig_heterogeneity']
)

#Rename columns
sex_pval_DF.columns = [f"{val}_{sex}" for val, sex in sex_pval_DF.columns]

# Reset index to turn geneID and cellType back into columns
sex_pval_DF = sex_pval_DF.reset_index()

#Rename columns for merging
sex_pval_DF = sex_pval_DF.rename(columns={'geneID':'gene_id'})

##Merge ES_diff results into metaanalysis summary files
for cell in ['CD4','CD8']:
    if cell == 'CD4':
        main_df = CD4_DF
    elif cell == 'CD8':
        main_df = CD8_DF
    
    #Merge ES_diff results
    sex_pval_DF2 = sex_pval_DF[sex_pval_DF['cellType'] == cell]
    
    merge_DF= pd.merge(main_df,sex_pval_DF2,on=['gene_id','cellType'], how = 'outer', indicator = 'bysex_merge_flag')
    
    merge_DF.to_csv(os.path.join(ind,f"metaanalysis_summary_{cell}_w_ES_diff_w_bysex.csv"), index = False)
    
    counts = merge_DF['bysex_merge_flag'].value_counts()
    print(counts)