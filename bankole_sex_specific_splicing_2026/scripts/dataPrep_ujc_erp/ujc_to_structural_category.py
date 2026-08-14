#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 13:51:41 2025

@author: nkeil
"""

import pandas as pd
import os


# Set paths to directories
rmg_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables"

#Set paths to each species ujc files
mel_file = os.path.join(rmg_dir,"rmg_mel_comb_TRs_facet_sex_ujc_counts.csv")
sim_file = os.path.join(rmg_dir,"rmg_sim_comb_TRs_facet_sex_ujc_counts.csv")
yak_file = os.path.join(rlr_dir,"TO_dyak_comb_TRs_facet_sex_ujc_counts.csv")
san_file = os.path.join(rlr_dir,"TO_dsan_comb_TRs_facet_sex_ujc_counts.csv")
ser_file = os.path.join(axk_dir,"kopp_dser_comb_TRs_facet_sex_ujc_counts.csv")

#Make dictionary linking species to filepath
species_files = {
    "mel": mel_file,
    "sim": sim_file,
    "yak": yak_file,
    "san": san_file,
    "ser": ser_file
}

for species, file_path in species_files.items():
    # Read the CSV file
    
    ujc_df = pd.read_csv(file_path)
    
    #Select only ujcs in annotated genes
    #!!! WE DECIDED we want structural categories for all UJCs
    #ujc_df = ujc_df[ujc_df['flag_annotated_gene'] == 1]
    
    ujc_cat_df = ujc_df[['jxnHash', 'structural_category']].drop_duplicates()
    
    cnt_df = ujc_cat_df.groupby('jxnHash').size().reset_index(name='category_count')
    
    #Merge category count into ujc_cat_df
    ujc_cat_df = pd.merge(ujc_cat_df, cnt_df, how='outer', on='jxnHash')
    
    # Set structural_category to "ambiguous" where category_count > 1
    ujc_cat_df.loc[ujc_cat_df['category_count'] > 1, 'structural_category'] = "ambiguous"
    
    #Drop the category_count column
    ujc_cat_df = ujc_cat_df.drop(columns=['category_count'])
    
    #Drop duplicates again
    ujc_cat_df = ujc_cat_df[['jxnHash', 'structural_category']].drop_duplicates()
    
    # Set output path
    output_file = os.path.join(outd, f"{species}_ujc_to_structural_category.csv")

    # Save to CSV with index=False to avoid adding an extra index column
    ujc_cat_df.to_csv(output_file, index=False)

    # # Counting the number of unique jxnHashes with >1 category count
    # num_gt_1 = cnt_df[cnt_df['category_count'] > 1]['jxnHash'].nunique()

    # # Counting the number of unique jxnHashes with =1 category count
    # num_eq_1 = cnt_df[cnt_df['category_count'] == 1]['jxnHash'].nunique()

    # # Dataframe with counts
    # result_df = pd.DataFrame({
    #     'Category Count Condition': ['>1', '=1'],
    #     'Unique jxnHash Count': [num_gt_1, num_eq_1]
    #     })
    
    




