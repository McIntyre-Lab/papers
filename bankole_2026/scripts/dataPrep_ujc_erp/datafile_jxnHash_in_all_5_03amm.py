#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:33:36 2025

@author: ammorse
"""

import pandas as pd
import gc

# Define file paths and species variable
base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

# Import the dmel6
data_dmel6 = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_dmel6_w_annoFlag_ERP_ESP_info.csv")
data2_dmel6 = data_dmel6[data_dmel6['flag_jxnHash_in_fiveSpecies_full_anno'] == 1]
print(data2_dmel6.columns)

## find variable for 1:1 species to dmel6_jxnHash

# rename columns
prefix = "dmel6_"
data2_dmel6 = data2_dmel6.rename(columns={
    col: f"{prefix}{col}" if col not in ['cat_dmel6_transcriptID'] else col
    for col in data2_dmel6.columns   
})
print(data2_dmel6.columns)


## (1) UJC in mel, sim and ser
species_list = ['dsim2', 'dser1']
dataframes = [data2_dmel6]

for species in species_list:
    # Import the comparison species data
    data = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info.csv")
    data2 = data[data['flag_jxnHash_in_fiveSpecies_full_anno'] == 1]   
    print(data2.columns)
    data2 = data2.rename(columns={'dmel6_jxnhash': 'dmel6_jxnHash'})
    data2 = data2.drop(['dmel6_geneID', 'cat_dmel6_transcriptID'], axis=1)
    
    # add prefix
    prefix = f"{species}_"
    data2 = data2.rename(columns={
        col: f"{prefix}{col}" if col not in ['dmel6_jxnHash', 'dmel6_geneID', 'cat_dmel6_transcriptID'] else col
        for col in data2.columns   
    })    
    print(data2.columns)
    # Append the processed data to the dataframes list
    dataframes.append(data2)
    
# Merge all dataframes on 'dmel6_jxnHash'
merged_df = dataframes[0]  # Start with the first dataframe (dmel6)
for i, df in enumerate(dataframes):
    print(f"DataFrame {i} columns:")
    print(df.columns.tolist())

for i, df in enumerate(dataframes[1:], start=1):  # Start from the second dataframe
    species = species_list[i-1]  # Get the species name for the current dataframe
    merged_df = pd.merge(merged_df, df, on='dmel6_jxnHash', how="inner")
    
# Check the merged dataframe
print(merged_df.head())
num_rows = merged_df.shape[0]
print(f"num rows after merging and filtering:  {num_rows}")
    # 9531 jxnHash in all 5 species.
num_uniq_geneID = merged_df['dmel6_geneID'].nunique()
print(f"num genes after merging and filtering: {num_uniq_geneID}")    

# include f and m counts for yak and san in dataset
data_dsan = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_dsan1_w_annoFlag_ERP_ESP_info.csv")
data_dsan = data_dsan.rename(columns={'dmel6_jxnhash': 'dmel6_jxnHash'})
print(data_dsan.columns)
keepers = ['jxnHash', 'geneID', 'rawCnts_dsan_F_rep1', 'rawCnts_dsan_M_rep1', 'dmel6_jxnHash']
data_dsan = data_dsan[keepers]
prefix = "dsan1_"
data_dsan = data_dsan.rename(columns={
     col: f"{prefix}{col}" if col not in ['dmel6_jxnHash'] else col
     for col in data_dsan.columns   
})    
print(data_dsan.columns)

data_dyak = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_dyak2_w_annoFlag_ERP_ESP_info.csv")
data_dyak = data_dyak.rename(columns={'dmel6_jxnhash': 'dmel6_jxnHash'})
print(data_dyak.columns)
keepers = ['jxnHash', 'geneID', 'rawCnts_dyak_F_rep1', 'rawCnts_dyak_M_rep1', 'dmel6_jxnHash']
data_dyak = data_dyak[keepers]
prefix = "dyak2_"
data_dyak = data_dyak.rename(columns={
     col: f"{prefix}{col}" if col not in ['dmel6_jxnHash'] else col
     for col in data_dyak.columns   
})    
print(data_dyak.columns)

# merge in dsan and dyak data
close = pd.merge(
    merged_df, 
    data_dsan, 
    on=["dmel6_jxnHash"], 
    how="outer", 
    indicator='merge_check')
print(close['merge_check'].value_counts(dropna=False).sort_index())
    #merge_check
    #left_only       1424
    #right_only    601961
    #both            9436
closer = close[close['merge_check'].str.contains('left_only|both', na=False)]

almost = pd.merge(
    closer, 
    data_dyak, 
    on=["dmel6_jxnHash"], 
    how="outer", 
    indicator='merge_check2')
print(almost['merge_check2'].value_counts(dropna=False).sort_index())
    #merge_check2
    #left_only       1745
    #right_only    512146
    #both            9950
almoster = almost[almost['merge_check'].str.contains('left_only|both', na=False)]
almoster = almoster.drop(['merge_check', 'merge_check2'], axis=1)

# save as csv 
almoster.to_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim_dser.csv", index=False)    

# clean up 
del data, data2, merged_df, dataframes, close, closer, almost, almoster
gc.collect()

## (2) UJC in mel and sim
species_list = ['dsim2']
dataframes = [data2_dmel6]

for species in species_list:
    # Import the comparison species data
    data = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info.csv")
    data2 = data[data['flag_jxnHash_in_fiveSpecies_full_anno'] == 1]   
    print(data2.columns)
    data2 = data2.rename(columns={'dmel6_jxnhash': 'dmel6_jxnHash'})
    data2 = data2.drop(['dmel6_geneID', 'cat_dmel6_transcriptID'], axis=1)
    
    # add prefix
    prefix = f"{species}_"
    data2 = data2.rename(columns={
        col: f"{prefix}{col}" if col not in ['dmel6_jxnHash', 'dmel6_geneID', 'cat_dmel6_transcriptID'] else col
        for col in data2.columns   
    })    
    print(data2.columns)
    # Append the processed data to the dataframes list
    dataframes.append(data2)
    
# Merge all dataframes on 'dmel6_jxnHash'
merged_df = dataframes[0]  # Start with the first dataframe (dmel6)

for i, df in enumerate(dataframes[1:], start=1):  # Start from the second dataframe
    species = species_list[i-1]  # Get the species name for the current dataframe
    merged_df = pd.merge(merged_df, df, on='dmel6_jxnHash', how="inner")
    
# Check the merged dataframe
print(merged_df.head())
num_rows = merged_df.shape[0]
print(f"num rows after merging and filtering:  {num_rows}")

num_uniq_geneID = merged_df['dmel6_geneID'].nunique()
print(f"num genes after merging and filtering: {num_uniq_geneID}")    

# merge in dsan and dyak data
close = pd.merge(
    merged_df, 
    data_dsan, 
    on=["dmel6_jxnHash"], 
    how="outer", 
    indicator='merge_check')
print(close['merge_check'].value_counts(dropna=False).sort_index())
    #merge_check
    #left_only       4214
    #right_only    596555
    #both            13870
closer = close[close['merge_check'].str.contains('left_only|both', na=False)]

almost = pd.merge(
    closer, 
    data_dyak, 
    on=["dmel6_jxnHash"], 
    how="outer", 
    indicator='merge_check2')
print(almost['merge_check2'].value_counts(dropna=False).sort_index())
    #merge_check2
    #left_only       4815
    #right_only    507061
    #both            14047
almoster = almost[almost['merge_check'].str.contains('left_only|both', na=False)]
almoster = almoster.drop(['merge_check', 'merge_check2'], axis=1)

# save as csv 
almoster.to_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim.csv", index=False)    
   

# 3) all 5 species 
# clean up befor next
del data, data2, merged_df, dataframes
gc.collect()

species_list = ['dsim2', 'dser1', 'dsan1', 'dyak2']
dataframes = [data2_dmel6]

for species in species_list:
    # Import the comparison species data
    data = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info.csv")
    data2 = data[data['flag_jxnHash_in_fiveSpecies_full_anno'] == 1]   
    print(data2.columns)
    data2 = data2.rename(columns={'dmel6_jxnhash': 'dmel6_jxnHash'})
    data2 = data2.drop(['dmel6_geneID', 'cat_dmel6_transcriptID'], axis=1)
    
    # add prefix
    prefix = f"{species}_"
    data2 = data2.rename(columns={
        col: f"{prefix}{col}" if col not in ['dmel6_jxnHash' ] else col
        for col in data2.columns   
    })    
    print(data2.columns)
    # Append the processed data to the dataframes list
    dataframes.append(data2)
    
# Merge all dataframes on 'dmel6_jxnHash'
merged_df = dataframes[0]  # Start with the first dataframe (dmel6)

for i, df in enumerate(dataframes[1:], start=1):  # Start from the second dataframe
    species = species_list[i-1]  # Get the species name for the current dataframe
    merged_df = pd.merge(merged_df, df, on='dmel6_jxnHash', how="inner")
    
# Check the merged dataframe
print(merged_df.head())
num_rows = merged_df.shape[0]
print(f"num rows after merging and filtering:  {num_rows}")
    # 9530 jxnHash in all 5 species.
num_uniq_geneID = merged_df['dmel6_geneID'].nunique()
print(f"num genes after merging and filtering: {num_uniq_geneID}")    
    # 4982 genes
    
# save as csv 
merged_df.to_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim_dser_dsan_dyak.csv", index=False)    




## some counting 
dmel_dsim = pd.read_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim.csv")
dmel_dsim_dser = pd.read_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim_dser.csv")
all5 = pd.read_csv(f"{base_path}/Tables/jxnHash_in_dmel_dsim_dser_dsan_dyak.csv")

# mel and sim counts
num_rows = dmel_dsim.shape[0]
print(f"num jxnhash in both dmel dsim:  {num_rows}")
    # 18862 jxnHash.
num_uniq_geneID = dmel_dsim['dmel6_geneID'].nunique()
print(f"num dmel6_geneID in mel and sim: {num_uniq_geneID}")    
    # 8679

# mel sim ser    
num_rows = dmel_dsim_dser.shape[0]
print(f"num jxnhash in both dmel dsim dser:  {num_rows}")
    # 11695 jxnHash.
num_uniq_geneID = dmel_dsim_dser['dmel6_geneID'].nunique()
print(f"num dmel6_geneID in mel and sim ser: {num_uniq_geneID}")    
    # 5756
    
# all 5    
num_rows = all5.shape[0]
print(f"num jxnhash in all 5 :  {num_rows}")
    # 13658 jxnHash.
num_uniq_geneID = all5['dmel6_geneID'].nunique()
print(f"num dmel6_geneID in all 5: {num_uniq_geneID}")    
    # 7800


