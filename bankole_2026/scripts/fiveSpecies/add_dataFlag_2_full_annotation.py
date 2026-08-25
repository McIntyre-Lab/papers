#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 10:07:12 2025

@author: ammorse

add flag to annotation files if jxnHash is in datafile

"""

import pandas as pd
import gc 
gc.collect()

# Define file paths and species variable
base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

species_list = ['dmel6', 'dsim2', 'dser1', 'dsan1', 'dyak2']
species_data = {}

for species in species_list:
    
    # import anno file
    anno = pd.read_csv(f"{base_path}/submission/supplementary/fiveSpecies_{species}_full_annotation.csv", low_memory=False)
    
    # import dmel datafile
    data = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv", low_memory=False)
    # create data flag, keep dmel6_jxnHash and flag
    data['flag_jxnHash_in_data'] = 1
    #print(data.columns)
    data = data.rename(columns={'jxnHash': f'{species}_jxnHash'})
    data = data[[f'{species}_jxnHash', 'flag_jxnHash_in_data']]
    data_flag = data.drop_duplicates(subset=[f'{species}_jxnHash', 'flag_jxnHash_in_data'])

    # merge data_flag into anno 
    anno_flag = pd.merge(
        data_flag,
        anno,
        on = f'{species}_jxnHash',
        how = 'outer',
        indicator = 'merge')
    print(f"{species} merge results:", anno_flag['merge'].value_counts(dropna=False).sort_index())
        #merge
        #left_only     2264256
        #right_only      53314
        #both            22672

    # keep if in right_only (annotation) and both (flag)
    anno_flag = anno_flag[anno_flag['merge'].isin(['right_only', 'both'])]
    # make nan 0s and convert to int
    anno_flag['flag_jxnHash_in_data'] = anno_flag['flag_jxnHash_in_data'].fillna(0).astype(int)
    print(f"{species} flag counts:", anno_flag['flag_jxnHash_in_data'].value_counts(dropna=False).sort_index())
    # Print total row count
    print(f"Total jxnHash in {species} annotation: {len(anno_flag)}")
    
    # Drop the 'merge' column
    anno_flag = anno_flag.drop(columns=['merge'])

    # Reorder columns to have f'{species}_jxnHash' and 'flag_jxnHash_in_data' first
    column_order = [f'{species}_jxnHash', 'flag_jxnHash_in_data'] + \
                   [col for col in anno.columns if col != f'{species}_jxnHash']

    # Reorder the dataframe columns
    anno_flag = anno_flag[column_order]
    
    ## save to csv
    anno_flag.to_csv(f"{base_path}/submission/supplementary/fiveSpecies_{species}_full_annotation_w_dataFlag.csv", index=False)    


"""
dmel6 flag counts: flag_jxnHash_in_data
0    53314
1    22672
Total jxnHash in dmel6 annotation: 75986

dsim2 flag counts: flag_jxnHash_in_data
0    53926
1    23355
Total jxnHash in dsim2 annotation: 77281

dser1 flag counts: flag_jxnHash_in_data
0    52583
1    19689
Total jxnHash in dser1 annotation: 72272

dsan1 flag counts: flag_jxnHash_in_data
0    55523
1    21500
Total jxnHash in dsan1 annotation: 77023

dyak2 flag counts: flag_jxnHash_in_data
0    59103
1    17937
Total jxnHash in dyak2 annotation: 77040
"""