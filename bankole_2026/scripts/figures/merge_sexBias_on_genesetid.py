#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 11:27:15 2026

@author: mgaran
"""
import pandas as pd

BASE_DIR    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY_DIR = f"{BASE_DIR}/zenodo/summary_files"
TABLES_DIR  = f"{BASE_DIR}/Tables"
OUTPUT_FILE = f"{TABLES_DIR}/gene_sexBias_merged_on_genesetid.csv"

sex_bias_species_list = ['dmel', 'dsim', 'dser']
analyzed_categories = ['F', 'M', 'B', 'unbiased']

merged_df = None

for sp in sex_bias_species_list:

    # Gene summary is already unique on geneID and already contains genesetid.
    df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_{sp}.csv",
        usecols=['geneID', 'genesetid', 'sexBias', 'num_ujc', 'num_ism_ujc'],
        low_memory=False,
    )

    print(f"\n{sp}")
    print(f"Gene summary rows: {len(df):,}")
    print(f"Unique geneIDs: {df['geneID'].nunique():,}")
    print(f"Unique genesetids: {df['genesetid'].nunique():,}")

    # Flag genes that are analyzable for sex bias.
    df['flag_analyzable'] = df['sexBias'].isin(analyzed_categories).astype(int)

    # B-biased genes need at least 2 non-ISM UJCs for Panel C.
    df['flag_B_multi'] = (
        (df['sexBias'] == 'B') &
        ((df['num_ujc'] - df['num_ism_ujc']) >= 2)
    ).astype(int)

    # Count genes of each sex-bias category within each genesetid.
    df['flag_M'] = (df['sexBias'] == 'M').astype(int)
    df['flag_F'] = (df['sexBias'] == 'F').astype(int)
    df['flag_B'] = (df['sexBias'] == 'B').astype(int)
    df['flag_unbiased'] = (df['sexBias'] == 'unbiased').astype(int)

    geneset_df = df.groupby('genesetid').agg(
        num_geneID=('geneID', 'nunique'),
        num_analyzable=('flag_analyzable', 'sum'),
        num_M=('flag_M', 'sum'),
        num_F=('flag_F', 'sum'),
        num_B=('flag_B', 'sum'),
        num_B_multi=('flag_B_multi', 'sum'),
        num_unbiased=('flag_unbiased', 'sum'),
    ).reset_index()

    geneset_df = geneset_df.rename(columns={
        'num_geneID': f'{sp}_num_geneID',
        'num_analyzable': f'{sp}_num_analyzable',
        'num_M': f'{sp}_num_M',
        'num_F': f'{sp}_num_F',
        'num_B': f'{sp}_num_B',
        'num_B_multi': f'{sp}_num_B_multi',
        'num_unbiased': f'{sp}_num_unbiased',
    })

    # Outer merge each species on genesetid.
    if merged_df is None:
        merged_df = geneset_df.copy()
    else:
        merged_df = pd.merge(
            merged_df,
            geneset_df,
            on='genesetid',
            how='outer',
            indicator='merge_check',
        )
        print("Merge check:")
        print(merged_df['merge_check'].value_counts(dropna=False))
        merged_df = merged_df.drop(columns='merge_check')

# A genesetid can be absent from one species. Treat its gene counts as zero.
count_cols = [c for c in merged_df.columns if c != 'genesetid']
merged_df[count_cols] = merged_df[count_cols].fillna(0).astype(int)

merged_df.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved {len(merged_df):,} genesetids to {OUTPUT_FILE}")
