#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 01:58:11 2025

@author: nkeil
"""

import pandas as pd
#import os

# Load your data
in_file = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/metaanalysis_summary_CD4_w_ES_diff_w_bysex.csv"
in_df = pd.read_csv(in_file, sep=",")  # Adjust separator if needed

# Identify all columns that start with 'flag_sig_'
flag_sig_cols = [col for col in in_df.columns if col.startswith("flag_sig_")]

# Initialize results
summary = []

# Total number of genes
total_genes = len(in_df)
total_immune_genes = in_df["flag_immune_gene"].sum()
total_t1d_genes = in_df["flag_t1d_gene"].sum()

for col in flag_sig_cols:
    sig_genes = in_df[in_df[col] == 1]
    total_count = len(sig_genes)
    immune_count = sig_genes['flag_immune_gene'].sum()
    t1d_count = sig_genes['flag_t1d_gene'].sum()
    
    
    
    summary.append({
    "flag_sig_column": col,
    "n_total_sig_genes": total_count,
    "prop_total_sig_genes": total_count / total_genes,
    "total_genes": total_genes,
    "n_immune_sig_genes": immune_count,
    "prop_immune_sig_genes": immune_count / total_immune_genes ,
    "total_immune_genes": total_immune_genes,
    "n_t1d_sig_genes": t1d_count,
    "prop_t1d_sig_genes": t1d_count / total_t1d_genes,
    "total_t1d_genes": total_t1d_genes
})

# Create a DataFrame for easy viewing
summary_df = pd.DataFrame(summary)

#Write csv
summary_df.to_csv(
    "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/CD4_metaanalysis_counts.csv",
    index=False
)