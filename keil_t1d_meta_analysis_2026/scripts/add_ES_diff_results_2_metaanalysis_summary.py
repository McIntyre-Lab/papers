#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 20 15:07:11 2025

@author: nkeil
"""

import os 
import pandas as pd

ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts"

CD4_DF=pd.read_csv(os.path.join(ind, "metaanalysis_summary_CD4.csv"))
CD8_DF=pd.read_csv(os.path.join(ind, "metaanalysis_summary_CD8.csv"))

es_diff_DF=pd.read_csv(os.path.join(ind, "meta_analysis_sex_diff_summary/meta_analysis_pvals_ES_diff_heterogenity.csv"))
#Rename columns for merging
es_diff_DF = es_diff_DF.rename(columns={'geneID':'gene_id','pval_heterogeneity': 'pval_het_sex_diff', 'flag_sig_heterogeneity': 'flag_sig_het_sex_diff'})

#Merge ES_diff results into metaanalysis summary files
for cell in ['CD4','CD8']:
    if cell == 'CD4':
        main_df = CD4_DF
    elif cell == 'CD8':
        main_df = CD8_DF
    
    #Merge ES_diff results
    es_diff_DF2 = es_diff_DF[es_diff_DF['cellType'] == cell]
    
    merge_DF= pd.merge(main_df,es_diff_DF2,on=['gene_id','cellType'], how = 'outer', indicator = "ES_diff_merge_ind")
    
    merge_DF.to_csv(os.path.join(ind,f"metaanalysis_summary_{cell}_w_ES_diff.csv"), index=False)
    
    counts = merge_DF["ES_diff_merge_ind"].value_counts()
    print(counts)
    
    
    