#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 14:58:45 2025

@author: ammorse

Add netanya’s structural category information to each:
datafile_jxnHash_{species}_w_ERP_ESP_info.csv

"""

import pandas as pd


base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

# Define file paths and species variable
species_dct = {"dmel6":"mel", "dsim2":"sim", "dser1":"ser", "dyak2":"yak", "dsan1":"san"}
#species_dct = {"dyak2":"yak"}

for species,species_data in species_dct.items():
    print(f"processing species: {species}")
   
    # import datafile                            
    datafile_jxnHash_species = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/old2/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_02amm.csv")
    print(datafile_jxnHash_species.columns)
    
    # import structural category info
    category = pd.read_csv(f"{base_path}/Tables/{species_data}_ujc_to_structural_category.csv")   
    print(category.columns)
    
    data_w_cat = pd.merge(
        datafile_jxnHash_species,
        category,
        on = 'jxnHash',
        how = 'outer',
        indicator = 'merge')
    print(f"{species} merge results:", data_w_cat['merge'].value_counts(dropna=False).sort_index())
        #dyak2 merge results: merge
        #left_only      11338
        #right_only         0
        #both          508340
    # keep if left only and both
    data_w_cat = data_w_cat[data_w_cat['merge'].isin(['left_only', 'both'])]
    # Drop the 'merge' column
    data_w_cat = data_w_cat.drop(columns=['merge'])
    print(data_w_cat.columns)
    print(f"{species} category counts:", data_w_cat['structural_category'].value_counts(dropna=False).sort_index())   
    
    # Reorder columns to have f'{species}_jxnHash' and 'flag_jxnHash_in_data' first
    first_cols = ['geneID','ERP', 'jxnHash', 'structural_category'] 
    other_cols = [col for col in data_w_cat.columns if col not in first_cols]

    # Reorder the dataframe columns
    data_w_cat = data_w_cat[first_cols + other_cols]
    
    ## save to csv
    data_w_cat.to_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_strCat.csv", index=False)    

## check no missing struc cat
    missing_values = data_w_cat['structural_category'].isnull().sum()

# Print the result
    print(f"Missing values in {species} 'structural_category' column: {missing_values}")