#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 16:15:32 2025

@author: nkeil
"""

import pandas as pd
import numpy as np
import os

#Set paths
put_transcript_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list_w_ERP.csv'
mel_ERP_info_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/fiveSpecies_2_dmel6_anno_files/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_infoERP.csv'

#Load transcript data
transcript_df = pd.read_csv(put_transcript_file, index_col=0)

#Dictionary linking species and tag
species_to_tag = {
    'mel': 'dmel6',
    'sim': 'dsim2',
    'yak': 'dyak2',
    'san': 'dsan1',
    'ser': 'dser1'
}

def add_flag_ERP_in_mel_anno(transcript_df, mel_erp_file):
    """
    Adds a column 'flag_ERP_in_mel_anno' to transcript_df that flags whether the ERP binary pattern
    (sign-stripped) exists for the corresponding mel_geneID in the mel ERP annotation file.

    Parameters:
    - transcript_df (pd.DataFrame): Must contain 'ERP' and 'mel_geneID' columns
    - mel_erp_file (str): Path to the mel ERP annotation CSV

    Returns:
    - pd.DataFrame: Copy of input with new column 'flag_ERP_in_mel_anno' (0 or 1)
    """
    transcript_df = transcript_df.copy()

    # Load mel ERP annotation and normalize ERP (strip +_ / -_)
    mel_erp_df = pd.read_csv(mel_erp_file, usecols=['jxnHash', 'geneID', 'ERP'])
    mel_erp_df['ERP_binary'] = mel_erp_df['ERP'].str.extract(r'[_-](\d+)$')  # pull only binary pattern

    # Build mapping: geneID → set of ERP binary patterns
    gene_to_erp_patterns = (
        mel_erp_df.groupby('geneID')['ERP_binary']
        .apply(set)
        .to_dict()
    )

    # Normalize ERP in transcript_df too (strip +_ or -_)
    transcript_df['ERP_binary'] = transcript_df['ERP'].str.extract(r'[_-](\d+)$')

    # Define flag: does the pattern exist for this gene?
    def check_erp_match(row):
        gene = row['mel_geneID']
        erp = row['ERP_binary']
        if pd.isna(gene) or pd.isna(erp):
            return 0
        return int(erp in gene_to_erp_patterns.get(gene, set()))

    transcript_df['flag_ERP_in_mel_anno'] = transcript_df.apply(check_erp_match, axis=1)

    # Drop the helper column
    transcript_df.drop(columns=['ERP_binary'], inplace=True)

    return transcript_df

def add_flag_ERP_in_anno_per_species(transcript_df, species_to_tag, base_erp_path):
    """
    Adds flag columns to transcript_df for each non-'mel' species indicating whether the ERP binary pattern
    for that row is present in that species' ERP annotation file for the corresponding gene (mel_geneID).

    Parameters:
    - transcript_df (pd.DataFrame): Must contain 'ERP', 'species', and 'mel_geneID' columns.
    - species_to_tag (dict): Maps species (e.g., 'sim') to tag (e.g., 'dsim2')
    - base_erp_path (str): Directory containing species-specific ERP files named like:
        'fiveSpecies_2_dmel6_ujc_er_vs_{tag}_ujc_infoERP.csv'

    Returns:
    - pd.DataFrame: Copy of transcript_df with new 'flag_ERP_in_{species}_anno' columns added.
    """
    transcript_df = transcript_df.copy()

    for species, tag in species_to_tag.items():
        if species == 'mel':
            continue  # Skip 'mel'

        # Path to species-specific ERP annotation file
        erp_file = os.path.join(
            base_erp_path,
            f'fiveSpecies_2_{tag}_anno_files/fiveSpecies_2_{tag}_ujc_er_vs_fiveSpecies_2_{tag}_ujc_infoERP.csv')

        try:
            erp_df = pd.read_csv(erp_file, usecols=['geneID', 'ERP'])
        except FileNotFoundError:
            print(f"⚠️ File not found for species '{species}': {erp_file}")
            transcript_df[f'flag_ERP_in_{species}_anno'] = np.nan
            continue

        # Strip sign from ERP to get binary pattern
        erp_df['ERP_binary'] = erp_df['ERP'].str.extract(r'[_-](\d+)$')
        gene_to_patterns = erp_df.groupby('geneID')['ERP_binary'].apply(set).to_dict()

        # Filter rows for this species
        species_mask = transcript_df['species'] == species
        working_df = transcript_df.loc[species_mask].copy()
        working_df['ERP_binary'] = working_df['ERP'].str.extract(r'[_-](\d+)$')

        def check_match(row):
            gene = row['associated_gene']
            erp = row['ERP_binary']
            if pd.isna(gene) or pd.isna(erp):
                return 0
            return int(erp in gene_to_patterns.get(gene, set()))

        # Apply and write back to main df
        transcript_df.loc[species_mask, f'flag_ERP_in_{species}_anno'] = working_df.apply(check_match, axis=1)

    return transcript_df

transcript_df = add_flag_ERP_in_mel_anno(transcript_df, mel_ERP_info_file)

erp_info_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations'

transcript_df = add_flag_ERP_in_anno_per_species(transcript_df, species_to_tag, erp_info_path)

transcript_df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list_w_ERP_flag_ERPinAnno.csv')