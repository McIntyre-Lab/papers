#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:53:39 2024

@author: ammorse
"""

"""
merge full annotation to datafile jxnhash 
"""

import pandas as pd

# Define file paths and species variable
species_dct = {"dmel6":"dmel", "dsim2":"dsim", "dser1":"dser", "dyak2":"dyak", "dsan1":"dsan"}
#species_dct = {"dsim2":"dsim"}
base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

for species,species_data in species_dct.items():
    print(f"processing species: {species}")
    # Import the species annotation CSV and datafile jxnHash
    species_anno = pd.read_csv(f"{base_path}/submission/supplementary/fiveSpecies_{species}_full_annotation_w_dataFlag.csv")
    species_anno = species_anno.rename(columns={f'{species}_jxnHash': 'jxnHash'})
    #print(species_anno.columns)                                            
    datafile_jxnHash_species = pd.read_csv(f"{base_path}/Tables/datafile_jxnHash_{species}.csv")
    #print(datafile_jxnHash_species.columns)

    # Import the infoERP and info ESP files - to merge to datafile
    species_ERP = pd.read_csv(f"{base_path}/erp_and_esp_output/fiveSpecies_2_{species}_ujc_er_vs_{species_data}_data_2_{species}_ujc_noMultiGene_infoERP.csv")
    #print(species_ERP.columns)
    species_ERP = species_ERP.rename(columns={
        'numExon': 'numExon_ERP',
        'flagDataOnlyExon': 'flagDataOnlyExon_ERP',
        'numDataOnlyExon':'numDataOnlyExon_ERP',
        'flagIR':'flagIR_ERP',
        'numIREvent': 'numIREvent_ERP',
        'flagReverseIR':'flagReverseIR_ERP',
        'flagNoSkip':'flagNoSkip_ERP',
        'flagNovel':'flagNovel_ERP',
        'flag5pFragment':'flag5pFragment_ERP',
        'flag3pFragment':'flag3pFragment_ERP',
        'flagIntrnlFrgmnt':'flagIntrnlFrgmnt_ERP',
        })
    species_ESP = pd.read_csv(f"{base_path}/erp_and_esp_output/fiveSpecies_2_{species}_ujc_es_vs_{species_data}_data_2_{species}_ujc_noMultiGene_infoESP.csv")
    #print(species_ESP.columns)
    species_ESP = species_ESP.rename(columns={
        'numExon': 'numExon_ESP',
        'flagDataOnlyExon': 'flagDataOnlyExon_ESP',
        'numDataOnlyExon': 'numDataOnlyExon_ESP'
        })
    
    # Select cols to keep from anno
    columns_to_keep = ['jxnHash', 'geneID', 'cat_dmel6_transcriptID', 
                       'flag_dsim202_2_dsim2_ujc',
                       'flag_dsimWXD_2_dsim2_ujc',
                       'flag_dyak21_2_dyak2_ujc',
                       'flag_dser11_2_dser1_ujc',
                       'flag_dmel650_2_dmel6_ujc',
                       'flag_dsan11_2_dsan1_ujc'
                       ]  
    species_anno = species_anno[columns_to_keep]
    species_anno_uniq = species_anno.drop_duplicates(subset=['jxnHash', 'geneID'])
    rows_before = len(species_anno)
    print(f"Rows {species} before: {rows_before}")  
    rows_after = len(species_anno_uniq)
    print(f"Rows {species} after: {rows_after}")
    
    # Merge the datasets
    merged = pd.merge(
        datafile_jxnHash_species, 
        species_anno, 
        on=["jxnHash", "geneID"], 
        how="outer", 
        indicator='merge_check')
    # print cols with _x or _y --> need rename if present
    cols_with_x_y = [col for col in merged.columns if col.endswith('_x') or col.endswith('_y')]
    print(cols_with_x_y)  # all good

    print(f"merge check for {species} jxnHash to full anno")
    print(merged['merge_check'].value_counts(dropna=False).sort_index())
        #merge_check
        #left_only     2264263  (data ujc only)
        #right_only      53321 (anno only )
        #both            22665 (data ujc and anno)
    #print(merged.columns)    
    # Drop rows if only in annotation (right_only)
    merged = merged[merged['merge_check'] != 'right_only']
    # create flag if data in annotation
    merged['flag_jxnHash_in_fiveSpecies_full_anno'] = merged['merge_check'].isin(['both']).astype(int)
    print(f"{species}")
    print(merged['flag_jxnHash_in_fiveSpecies_full_anno'].value_counts())
        #flag_in_fiveSpecies_full_anno
        #0          2264263
        #1          22665
    #print(merged.columns)
    # drop merge check col
    merged = merged.drop(columns=['merge_check'])

    ## merge in ERP and ESP info files
    merged2 = pd.merge(
        merged, 
        species_ERP, 
        on=["jxnHash", "geneID"], 
        how="outer", 
        indicator='merge_check')   
    print(merged2['merge_check'].value_counts())
        #merge_check
        #both          2286928
        #left_only           0
        #right_only          0
    # print cols with _x or _y --> need rename if present
    cols_with_x_y = [col for col in merged.columns if col.endswith('_x') or col.endswith('_y')]
    print(cols_with_x_y)  # all good  
    
    merged3 = pd.merge(
        merged2, 
        species_ESP, 
        on=["jxnHash", "geneID"], 
        how="outer", 
        indicator='merge_check2')   
    #print(merged3['merge_check2'].value_counts())
        #merge_check2
        #both          2286928
        #left_only           0
        #right_only          0
    cols_with_x_y = [col for col in merged.columns if col.endswith('_x') or col.endswith('_y')]
    print(cols_with_x_y)  # all good      
    
    # count number of rows before and after merging in fiveSpecies anno
    rows_before = len(datafile_jxnHash_species)
    print(f"Rows {species} jxnHash before merging: {rows_before}")    
    rows_after = len(merged3)
    print(f"Rows {species} after merging: {rows_after}")

    #print(merged3.columns)
    # drop merge indicator variables
    merged3 = merged3.drop(['merge_check', 'merge_check2'], axis=1)
    # check have all columns wanted
    #print(merged3.columns)
    
    # save to csv
    merged3.to_csv(f"{base_path}/submission/supplementary/datafiles/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info.csv", index=False)
    print(f"successfully processed species: {species}")
