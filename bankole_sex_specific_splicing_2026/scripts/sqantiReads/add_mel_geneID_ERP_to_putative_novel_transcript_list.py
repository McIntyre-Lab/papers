#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 15:03:00 2025

@author: nkeil
"""
import pandas as pd
import os
import numpy as np


# === Load inputs ===
put_transcript_list = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list.csv'
ortholog_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_5species_list_02ndk.csv'

# Load data
transcript_df = pd.read_csv(put_transcript_list, index_col=0)
ortholog_df = pd.read_csv(ortholog_file)

species_to_tag = {
    'mel': 'dmel6',
    'sim': 'dsim2',
    'yak': 'dyak2',
    'san': 'dsan1',
    'ser': 'dser1'
}

def add_mel_geneID(transcript_df, ortholog_df):
    """
    Add a column 'mel_geneID' to transcript_df using ortholog_df to map species-specific gene IDs to D. melanogaster gene IDs.
    
    Parameters:
    - transcript_df (pd.DataFrame): DataFrame containing 'species' and 'associated_gene' columns
    - ortholog_df (pd.DataFrame): DataFrame containing columns for gene IDs across species including 'dmel_geneID'
    
    Returns:
    - pd.DataFrame: A copy of transcript_df with an added 'mel_geneID' column (with np.nan for missing).
    """
    transcript_df = transcript_df.copy()
    transcript_df['mel_geneID'] = np.nan  # Set default to NaN

    for idx, row in transcript_df.iterrows():
        species = row['species']
        associated_gene = row['associated_gene']

        if species == 'mel':
            transcript_df.at[idx, 'mel_geneID'] = associated_gene
        else:
            ortholog_col = f"d{species}_geneID"
            if ortholog_col in ortholog_df.columns:
                match = ortholog_df[ortholog_df[ortholog_col] == associated_gene]
                if not match.empty:
                    transcript_df.at[idx, 'mel_geneID'] = match['dmel_geneID'].values[0]

    # Report number of missing mel_geneIDs by species
    missing_counts = transcript_df[transcript_df['mel_geneID'].isna()]['species'].value_counts()
    print("Number of transcripts without mel_geneID per species:")
    print(missing_counts)

    return transcript_df


def add_ERP(transcript_df, species_to_tag, base_erp_path):
    """
    Adds a single 'ERP' column to transcript_df by mapping species-specific jxnHash values using their ERP files.

    Parameters:
    - transcript_df (pd.DataFrame): Must contain 'jxnHash', 'species', and 'tag' columns.
    - species_to_tag (dict): Maps species keys (e.g., 'mel') to tag suffixes (e.g., 'dmel6').
    - base_erp_path (str): Directory path containing ERP files named like:
        'datafile_jxnHash_{tag}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv'

    Returns:
    - pd.DataFrame: Copy of input with a new 'ERP' column.
    """
    transcript_df = transcript_df.copy()
    transcript_df['ERP'] = None

    # Load all ERP mappings into a dict of dicts: {tag: {jxnHash: ERP}}
    erp_maps = {}
    for species, tag in species_to_tag.items():
        erp_file = os.path.join(
            base_erp_path,
            f'datafile_jxnHash_{tag}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv'
        )
        try:
            erp_df = pd.read_csv(erp_file, usecols=[f'{tag}_jxnHash', 'ERP'])
            erp_df = erp_df.rename(columns={f'{tag}_jxnHash': 'jxnHash'})
            erp_maps[tag] = erp_df.set_index('jxnHash')['ERP'].to_dict()
        except FileNotFoundError:
            print(f"⚠️ File not found: {erp_file}")
            erp_maps[tag] = {}

    # Apply ERP mapping row-wise using species-specific ERP map
    transcript_df['ERP'] = transcript_df.apply(
        lambda row: erp_maps.get(row['tag'], {}).get(row['jxnHash'], None),
        axis=1
    )

    return transcript_df

#Add mel geneID
#Note cannot link to mel gene ID for some jxnHashes
transcript_df = add_mel_geneID(transcript_df, ortholog_df)

#Add ERP
base_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/datafiles'
transcript_df = add_ERP(transcript_df, species_to_tag, base_path)

transcript_df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list_w_ERP.csv')




        
    
