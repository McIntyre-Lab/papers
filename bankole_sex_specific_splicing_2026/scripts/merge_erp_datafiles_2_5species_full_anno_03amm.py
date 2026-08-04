#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 10:25:17 2024

@author: ammorse
"""
"""
from 5species full - keep geneID and ERP_plus only, make uniq 
    create flag_in_anno = 1;
merge this erp_list into erp results 
    use flag_in_anno to split into +/- annotation
"""

import pandas as pd
# Define file paths and species variable

species_dct = {"dmel6":"dmel", "dsim2":"dsim", "dser1":"dser", "dyak2":"dyak", "dsan1":"dsan"}
#species_dct = {"dsim2":"dsim"}
 
base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

for species,species_data in species_dct.items():
    print(f"processing species: {species}")
    #species = "dmel6"
    #species_data = "dmel"
    
    # Import the species annotation CSV and data ERP file
    species_anno = pd.read_csv(f"{base_path}/submission/supplementary/fiveSpecies_{species}_full_annotation.csv")
    print(f"rows in {species} anno: {len(species_anno)}")
    datafile_erp_species = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_erp_{species}.csv", low_memory=False)
    print(f"rows in {species} datafile: {len(datafile_erp_species)}")
    
    # rename erp_plus col
    datafile_erp_species2 = datafile_erp_species.rename(columns={"erp_plus": "ERP_plus"})
    #print(datafile_erp_species2.columns)
    
    # in anno, keep geneID and erp_plus cols, make uniq and create flag
    anno_cols = ['geneID', 'ERP_plus']
    anno_keep = species_anno[anno_cols]
    anno_uniq = anno_keep.drop_duplicates().copy()
    print(f"rows in {species} anno_uniq: {len(anno_uniq)}")
    anno_uniq['flag_in_full_anno'] = 1
    
    # Merge anno_uniq into ERP results
    merged = pd.merge(
        datafile_erp_species2, 
        anno_uniq, 
        on=["geneID", "ERP_plus"], 
        how="outer", 
        indicator='merge_check')
    print(merged['merge_check'].value_counts(dropna=False).sort_index())
        #merge_check sim
        #left_only     725503  only in erp datafile
        #right_only     18154  only in anno file
        #both           29873  in both

    # Drop rows if only in annotation (right_only)
    merged = merged[merged['merge_check'] != 'right_only']          
 
    # if flag_in_full_anno is missing then make 0
    merged['flag_in_full_anno'] = merged['flag_in_full_anno'].fillna(0).astype(int)
    num_flag_1 = (merged['flag_in_full_anno'] == 1).sum()
    num_flag_0 = (merged['flag_in_full_anno'] == 0).sum()
    print(f"{species} in flag_in_full_anno = 1: {num_flag_1}")
    print(f"{species} not flag_in_full_anno = 0: {num_flag_0}")
        
    ## check to make sure var names are not overlapping
    cols_with_x_y = [col for col in merged.columns if col.endswith('_x') or col.endswith('_y')]
    print(f"for {species}: {cols_with_x_y}")  # all good 

    #print(merged.columns)
    merged = merged.drop(['merge_check'], axis=1)
    
    # save to csv
    merged.to_csv(f"{base_path}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag.csv", index=False)
    print(f"successfully processed species: {species}")
