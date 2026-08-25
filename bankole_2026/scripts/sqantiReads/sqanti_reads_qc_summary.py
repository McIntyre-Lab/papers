#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 23:19:54 2025

@author: nkeil
"""
import pandas as pd
import os
import re


# Set paths to directories
rmg_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"

def concat_csvs(*file_paths, axis=0, ignore_index=True):
    """
    Load multiple CSV files and concatenate them into a single DataFrame.
    
    Parameters:
    - *file_paths: Variable number of CSV file paths.
    - axis: The axis along which to concatenate. Default is 0 (stack rows).
            Use 1 to concatenate columns.
    - ignore_index: If True, the resulting DataFrame will have a new integer index.
    
    Returns:
    - concatenated_df: A DataFrame resulting from concatenating all loaded CSV files.
    """
    dataframes = [pd.read_csv(file) for file in file_paths]
    concatenated_df = pd.concat(dataframes, axis=axis, ignore_index=ignore_index)
    return concatenated_df

def extract_replicate(sample_id):
    """
    Extracts the number following 'rep' from a sample_id string.
    If no 'rep' is found, returns 1.
    """
    match = re.search(r'rep(\d+)', str(sample_id), re.IGNORECASE)
    return int(match.group(1)) if match else 1

#Concat length summary + category percentages
mel_file = os.path.join(rmg_dir,"rmg_mel_comb_TRs_facet_sex_length_summary_w_category_percs.csv")
sim_file = os.path.join(rmg_dir,"rmg_sim_comb_TRs_facet_sex_length_summary_w_category_percs.csv")
yak_file = os.path.join(rlr_dir,"TO_dyak_comb_TRs_facet_sex_length_summary_w_category_percs.csv")
san_file = os.path.join(rlr_dir,"TO_dsan_comb_TRs_facet_sex_length_summary_w_category_percs.csv")
ser_file = os.path.join(axk_dir,"kopp_dser_comb_TRs_facet_sex_length_summary_w_category_percs.csv")

length_df = concat_csvs(mel_file,sim_file,ser_file, yak_file,san_file)
#Drop sex column
length_df.drop(columns=['sex'], inplace=True)


#Concat artefact files
mel_file = os.path.join(rmg_dir,"rmg_mel_comb_TRs_facet_sex_err_counts.csv")
sim_file = os.path.join(rmg_dir,"rmg_sim_comb_TRs_facet_sex_err_counts.csv")
yak_file = os.path.join(rlr_dir,"TO_dyak_comb_TRs_facet_sex_err_counts.csv")
san_file = os.path.join(rlr_dir,"TO_dsan_comb_TRs_facet_sex_err_counts.csv")
ser_file = os.path.join(axk_dir,"kopp_dser_comb_TRs_facet_sex_err_counts.csv")

err_df = concat_csvs(mel_file,sim_file,ser_file, yak_file,san_file)

#Merge err_df and length_df
summary_df = pd.merge(err_df, length_df, on = "sampleID", how ="outer")

#Add a rep column
# Extract the last value from the 'sampleID' column
rep_value = summary_df['sampleID'].iloc[-1]

# Insert a new column 'rep' at position 1 (i.e. as the second column) with the constant rep_value
summary_df['replicate'] = summary_df['sampleID'].apply(extract_replicate)

# Reorder columns
cols = list(summary_df.columns)
cols.insert(1, cols.pop(cols.index('replicate')))
summary_df = summary_df[cols]

#Add perc_ prefix t structural category columns
summary_df.rename(columns=lambda col: f"perc_{col}" if col and col[0].isupper() else col, inplace=True)

summary_df['perc_FSM_ISM_NIC'] = summary_df['perc_FSM'] + summary_df['perc_ISM'] + summary_df['perc_NIC']

#Write summary_df to a csv file
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables"
summary_df.to_csv(os.path.join(outd,"sqanti_reads_qc_summary.csv"))