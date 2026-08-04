#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 15:38:44 2025

@author: ammorse

add flag_ERPanno to datafile_jxnHash

from 5species full anno, extract geneID and ERP columns, make uniq
merge geneID-ERP (L) to datafile jxnHash (R)
R-only = flag_ERPanno = 0
L-only = flag_ERPanno = 1

"""

import pandas as pd


base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

# Define file paths and species variable
species_dct = {"dmel6":"dmel", "dsim2":"dsim", "dser1":"dser", "dyak2":"dyak", "dsan1":"dsan"}
#species_dct = {"dyak2":"dyak"}

for species,species_data in species_dct.items():
    print(f"processing species: {species}")
    # Import the species annotation CSV and datafile jxnHash
    species_anno = pd.read_csv(f"{base_path}/submission/supplementary/fiveSpecies_{species}_full_annotation_w_dataFlag.csv")
    #print(species_anno.columns)     
    # keep geneID and ERP - make uniq
    species_anno2 = species_anno[['geneID', 'ERP']].drop_duplicates()
    
    # import datafile                            
    datafile_jxnHash_species = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info.csv")
    #print(datafile_jxnHash_species.columns)
    
    add_flag = pd.merge(
        species_anno2,
        datafile_jxnHash_species,
        on = ['geneID', 'ERP'],
        how = 'outer',
        indicator = 'merge')
    print(add_flag['merge'].value_counts(dropna=False).sort_index())
        #merge
        #left_only      22349
        #right_only    216884
        #both          302794  
    # Drop rows where _merge == 'left_only'
    add_flag = add_flag[add_flag['merge'] != 'left_only']

    # Create flag_anno where it's 1 if the merge result is 'both', otherwise 0
    add_flag['flag_ERPanno'] = add_flag['merge'].apply(lambda x: 1 if x == 'both' else 0)

    # Drop the '_merge' column for clarity
    add_flag = add_flag.drop(columns=['merge'])

    print(add_flag['flag_ERPanno'].value_counts(dropna=False).sort_index())
        #flag_ERPanno
        #0    216884
        #1    302794
        
    # save to CSV 
    add_flag.to_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_02amm.csv", index = False)

