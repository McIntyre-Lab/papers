#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 17:16:05 2025

@author: ammorse

count num exons in annotation and in data for each species

"""
import pandas as pd

base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

species_dct = {"dmel6":"dmel", "dsim2":"dsim", "dser1":"dser", "dsan1": "dsan", "dyak2": "dyak"}
#species_dct = {"dmel6":"dmel"}

for species,species_data in species_dct.items():
    print(f"processing species: {species}")
    
    #data_files = "/TB20/sex_specific_splicing/erp_and_esp_analysis"
    df_data = pd.read_csv(f"{base_path}/erp_output/fiveSpecies_2_{species}_ujc_er_vs_{species_data}_data_2_{species}_ujc_noMultiGene_infoERP.csv", low_memory=False)
    print(df_data.columns)
    keep = ['geneID', 'ERP', 'ERP_plus', 'numExon']
    df_data = df_data[keep].copy()    
    df_data = df_data.drop_duplicates()
    df_data['flag_ERP_1_exon'] = df_data['ERP'].str.fullmatch(r'[+-]_1').astype(int)
    df_data_mono = df_data[(df_data['flag_ERP_1_exon'] == 1) & (df_data['numExon'] == 1)]
   
    df_data_mono = df_data_mono.drop_duplicates()
    
    ## output as csv
    df_data_mono.to_csv(f"{base_path}/Tables/monoExon_{species}_in_data_infoERP.csv", index=False)   

    
    anno_files = f"{base_path}/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_{species}_anno_files"
    # mel exon num in anno
    df_anno = pd.read_csv(f"{anno_files}/fiveSpecies_2_{species}_ujc_er_vs_fiveSpecies_2_{species}_ujc_infoERP.csv", low_memory=False)
    #print(df_anno.columns)
    keep = ['geneID', 'ERP', 'ERP_plus', 'numExon']
    df_anno = df_anno[keep].copy()
    df_anno = df_anno.drop_duplicates()
    df_anno['flag_ERP_1_exon'] = df_anno['ERP'].str.fullmatch(r'[+-]_1').astype(int)
    df_anno_mono = df_anno[(df_anno['flag_ERP_1_exon'] == 1) & (df_anno['numExon'] == 1)]
    df_anno_mono = df_data_mono.drop_duplicates()

    ## output as csv
    df_anno_mono.to_csv(f"{base_path}/Tables/monoExon_{species}_in_anno_infoERP.csv", index=False)   

    
    