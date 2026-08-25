#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 13:10:50 2025

@author: ammorse

"""


import pandas as pd

base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

# import putative novel 
novel = pd.read_csv(f"{base_path}/sqanti_reads_analysis/putative_novel_transcript_list.csv", low_memory=False)
print(novel.columns)
# split novel by tag 
tags = novel.groupby('tag')
for species, group_df in tags:
    print(f"dataframe for tag = {species}")
    group_df['flag_putative_novel_jxnHash'] = 1  
    ## keep cols by name 
    keep = ['jxnHash', 'flag_putative_novel_jxnHash']
    group_df2 =group_df[keep] 
    group_df2 = group_df2.rename(columns={'jxnHash': f"{species}_jxnHash"}).reset_index(drop=True)
    print(group_df2.index)

    # import datafile for species
    df_species = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_strCat.csv", low_memory=False)
    #print(df_species.columns)
    df_species = df_species.rename(columns={'jxnHash': f"{species}_jxnHash"}).reset_index(drop=True)
    print(df_species.index)

    df_species_wFlag = pd.merge(
        df_species,
        group_df2,  # Merge using columns, not the index
        on=f"{species}_jxnHash",
        how='outer',
        indicator='merge'
    )

    print(df_species_wFlag['merge'].value_counts(dropna=False).sort_index())
        #left_only     2286875
        #right_only          0
        #both               53
    # keep if in left only and in both     
    df_species_wFlag2 = df_species_wFlag[df_species_wFlag['merge'].isin(['left_only', 'both'])]
    # if flag is missing then make 0 
    df_species_wFlag2['flag_putative_novel_jxnHash'] = df_species_wFlag2['flag_putative_novel_jxnHash'].fillna(0).astype(int)
    # Drop the 'merge' column
    df_species_wFlag2 = df_species_wFlag2.drop(columns=['merge'])
    
    ## save to csv
    df_species_wFlag2.to_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv", index=False)    

