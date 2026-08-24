#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 13:35:50 2026

@author: mgaran
"""

import pandas as pd

BASE     = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
OUTDIR = f"{BASE}/zenodo/summary_files"

SPECIES_LIST     = ['dmel', 'dsim', 'dser', 'dsan', 'dyak']
ANALYZED_SPECIES = ['dmel', 'dsim', 'dser']

G1 = 'F'
G2 = 'M'

def main():
    log_file = f"{BASE}/Tables/log_files/merge_component_summaries.log"

    with open(log_file, 'w') as f:
        f.write("Merging component summaries across all species\n\n")

    # add species prefix to columns from each individual species component summary
    all_components = set()
    species_dfs    = {}

    for species in SPECIES_LIST:
        file_path = f"{OUTDIR}/component_summary_{species}.csv"
        df = pd.read_csv(file_path, low_memory=False)
        all_components.update(df['componentID'])

        rename_dict = {col: f"{col}_{species}" for col in df.columns if col != 'componentID'}
        df = df.rename(columns=rename_dict)
        species_dfs[species] = df

    # turn series of components to dataframe before merging rows from each species component summary
    merged = pd.DataFrame({'componentID': sorted(all_components)})
    
    with open(log_file, 'a') as f:
        f.write("Merging component summaries across all species\n\n")
        f.write(f"Check that the table with only componentID columns is unique before any merges: {merged['componentID'].is_unique}")
            
    for species in SPECIES_LIST:
        merged_before = len(merged)
        merged = pd.merge(merged, species_dfs[species], on='componentID', how='outer', indicator='merge_check')

        # add num of components in all_components
        with open(log_file, 'a') as f:
            f.write(f"{species}:\n")
            f.write(f"  Rows before merge: {merged_before:,}\n")
            f.write(f"  Species rows:      {len(species_dfs[species]):,}\n")
            f.write("  Merge results:\n")
            f.write(f"{merged['merge_check'].value_counts(dropna=False).sort_index()}\n\n")

        merged = merged.drop('merge_check', axis=1)

    # Cross-species bias (analyzed species only)
    def get_pattern(row, col_suffix):
        pattern = ''
        for sp in ANALYZED_SPECIES:
            val = row.get(f'{col_suffix}_{sp}', '0')
            pattern += str(val) if pd.notna(val) else '0'
        return pattern

    merged['component_bias_pattern'] = merged.apply(
        lambda row: get_pattern(row, 'component_bias'), axis=1
    )

    INT_COLS = [
        'flag_simple', 'flag_monoExon',
        'num_UJC', 'num_ERPp', 'num_geneID',
        'num_UJC_w_data', 'num_ERPp_w_data',
        'num_geneID_w_data',
        f'num_UJC_bias_{G1}', f'num_UJC_bias_{G2}',
        'num_UJC_analyzable',
        'sumReads', f'sumReads_{G1}', f'sumReads_{G2}',
    ]

    for species in SPECIES_LIST:
        for col_base in INT_COLS:
            col = f'{col_base}_{species}'
            if col in merged.columns:
                merged[col] = merged[col].fillna(0).astype(int)

    output = f"{OUTDIR}/component_summary.csv"
    merged.to_csv(output, index=False)

if __name__ == '__main__':
    main()
