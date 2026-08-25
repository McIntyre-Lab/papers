#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 17:46:29 2025

@author: nkeil
"""
import pandas as pd

# Set up paths to ortholog files
dmel_2_dsim_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/useful_dmel_dsim_compare/dmel_orthologs_dsim_fb_2017_04.csv'
dmel_2_dyak_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_merge_dmel_dyak.csv'
dmel_2_dsan_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_merge_dmel_dsan.csv'
dmel_2_dser_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_merge_dmel_dser.csv'

# Load mel-dsim orthologs
dmel_2_dsim_df = pd.read_csv(dmel_2_dsim_file, usecols=['mel_geneID', 'sim_geneID'])

# Load dyak orthologs
dmel_2_dyak_df = pd.read_csv(dmel_2_dyak_file)
ortholog_df = pd.merge(dmel_2_dsim_df, dmel_2_dyak_df, left_on='mel_geneID', right_on='dmel_geneID', how='outer', indicator='yak_merge')
ortholog_df['mel_geneID'] = ortholog_df['mel_geneID'].fillna(ortholog_df['dmel_geneID'])
ortholog_df = ortholog_df.drop(columns=['ortholog_set', 'dmel_geneID'])

# Load dsan orthologs
dmel_2_dsan_df = pd.read_csv(dmel_2_dsan_file)
ortholog_df = pd.merge(ortholog_df, dmel_2_dsan_df, left_on='mel_geneID', right_on='dmel_geneID', how='outer', indicator='san_merge')
ortholog_df['mel_geneID'] = ortholog_df['mel_geneID'].fillna(ortholog_df['dmel_geneID'])
ortholog_df = ortholog_df.drop(columns=['ortholog_set', 'dmel_geneID'])

# Load dser orthologs
dmel_2_dser_df = pd.read_csv(dmel_2_dser_file)
ortholog_df = pd.merge(ortholog_df, dmel_2_dser_df, left_on='mel_geneID', right_on='dmel_geneID', how='outer', indicator='ser_merge')
ortholog_df['mel_geneID'] = ortholog_df['mel_geneID'].fillna(ortholog_df['dmel_geneID'])
ortholog_df = ortholog_df.drop(columns=['ortholog_set', 'dmel_geneID'])

# Format final DataFrame
ortholog_df = ortholog_df.rename(columns={
    'mel_geneID': 'dmel_geneID',
    'sim_geneID': 'dsim_geneID',
    'dyak_geneID': 'dyak_geneID',
    'dsan_geneID': 'dsan_geneID',
    'dser_geneID': 'dser_geneID'
}).copy()

# Add "LOC" prefix to numeric geneIDs
for col in ['dyak_geneID', 'dsan_geneID', 'dser_geneID']:
    mask = ortholog_df[col].notna()
    ortholog_df.loc[mask, col] = ortholog_df.loc[mask, col].astype(int).astype(str).apply(lambda x: f'LOC{x}')

# Ensure dmel and dsim IDs are strings
for col in ['dmel_geneID', 'dsim_geneID']:
    ortholog_df[col] = ortholog_df[col].astype(str)

# Save to file
ortholog_df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_5species_list_02ndk.csv', index=False)