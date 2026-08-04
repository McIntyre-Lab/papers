#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 11:45:18 2025

@author: ammorse
""" 
import pandas as pd
import numpy as np
   
# Define file paths and species variable

species_dct = {"dsan1":"dsan", "dyak2":"dyak"}
#species_dct = {"dsan1":"dsan"}
base_path = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

for species,species_data in species_dct.items():
    print(f"processing species: {species}")

    # Import the datafile for species
    df_species = pd.read_csv(f"{base_path}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag_flagNovel.csv")
    print(df_species.columns)

    # Convert to numeric, coerce errors to NaN
    df_species['flag_sex_bias_F'] = pd.to_numeric(df_species['flag_sex_bias_F'], errors='coerce')
    df_species['flag_sex_bias_M'] = pd.to_numeric(df_species['flag_sex_bias_M'], errors='coerce')

    df_species['newcol'] = df_species['flag_in_full_anno'].apply(
    lambda x: 'in_anno' if x == 1 else ('not_in_anno' if x == 0 else 'oops')
    )
    
    ## by gene, count sig ujc with F and M bias including category
    gene_summary = df_species.groupby(['geneID', 'newcol'])[
        ['flag_sex_bias_F', 'flag_sex_bias_M']].sum()
   
    # Pivot to wide format
    pivoted = gene_summary.unstack('newcol', fill_value=0)
    print(pivoted.columns)
    
    pivoted.columns = [
    f"num_erp_sexBias_{bias.replace('flag_bias_', '')}_{cat}"
    for bias, cat in pivoted.columns
    ]

    # Identify the columns for F and M
    f_cols = [col for col in pivoted.columns if '_F_' in col]
    m_cols = [col for col in pivoted.columns if '_M_' in col]
    
    # Check presence of at least one event across structural categories
    has_f = pivoted[f_cols].sum(axis=1) > 0
    has_m = pivoted[m_cols].sum(axis=1) > 0
    
    # Assign category
    pivoted['erp_sexBias'] = np.select(
        [
            has_f & ~has_m,   # F only
            ~has_f & has_m,   # M only
            has_f & has_m     # Both
        ],
        ['F', 'M', 'B'],
        default='.'
    )

    # geneID is index = reset it
    pivoted = pivoted.reset_index()
    # Get all current columns, define order and reorder
    all_cols = list(pivoted.columns)
    first_cols = ['geneID', 'erp_sexBias']    
    remaining_cols = [col for col in all_cols if col not in first_cols]
    pivoted = pivoted[first_cols + remaining_cols]
    print(pivoted.columns)
        
    for col in pivoted.columns:
        if col == 'geneID':
            continue  # Skip geneID column
        if pivoted[col].notna().any():
            print(f"\n=== Value counts for: {col} ===")
            print(pivoted[col].value_counts(dropna=False))
    print(pivoted['erp_sexBias'].value_counts(dropna=False))
        
    pivoted.to_csv(f"{base_path}/Tables/roz_stuff/geneCnts_num_bias_erp_by_Anno_{species}.csv", index=False)
