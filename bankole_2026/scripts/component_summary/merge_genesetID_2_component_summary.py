#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:55:21 2026

@author: mgaran
"""

import pandas as pd

PROJ   = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
LOG    = f"{PROJ}/Tables/log_files/merge_genesetID_to_component_summary.log"

NODES_FILE   = f"{PROJ}/zenodo/fiveSpecies_supporting_files/network_files/nodes_with_geneset.csv"
COMP_FILE    = f"{PROJ}/zenodo/summary_files/component_summary.csv"
OUTPUT_FILE  = f"{PROJ}/zenodo/summary_files/component_summary.csv"

def main():
    nodes_df = pd.read_csv(NODES_FILE, dtype=str)

    ## Keep only component nodes (nodeid starts with "c")
    comp_nodes = nodes_df[nodes_df['nodeid'].str.startswith('c')].copy()
    # remove 'c' at the beggining of the nodeid which was previously added to
    # the componentID variable to fix mixxed dtype error when creating a
    # network with both componentID (integers) and geneID (strings)
    comp_nodes['componentID'] = comp_nodes['nodeid'].str[1:].astype(int)
    comp_nodes = comp_nodes[['componentID', 'genesetid']].drop_duplicates()
    comp_nodes['genesetid'] = comp_nodes['genesetid'].astype(int)

    comp_df = pd.read_csv(COMP_FILE, low_memory=False)

    merged = pd.merge(comp_df, comp_nodes, on='componentID', how='outer', indicator='merge_check')
    with open(LOG, 'w') as log:
        log.write("Merge genesetID to merged_fiveSpecies_component_summary\n\n")
        log.write(f"Component summary rows:        {comp_df.shape[0]:,}\n")
        log.write(f"Component nodes with genesetid: {comp_nodes.shape[0]:,}\n")
        log.write(f"check component nodes unique on componentID: {comp_nodes['componentID'].is_unique}\n\n")
        log.write("Merge results (component_summary left join genesetid):\n")
        log.write(f"{merged['merge_check'].value_counts(dropna=False).sort_index()}\n\n\n")
        log.write(f"Count the total number of unique genesetid: {merged['genesetid'].nunique()}\n")
        log.write(f"Count the number of unique genesetid with data: {merged.loc[merged['genesetid'].notna(), 'genesetid'].nunique()}\n")

    merged = merged[merged['merge_check'].isin(['left_only', 'both'])]
    merged = merged.drop('merge_check', axis=1)

    # Move genesetid to second column (after componentID)
    cols = merged.columns.tolist()
    cols.remove('genesetid')
    cols.insert(1, 'genesetid')
    merged = merged[cols]
    merged.to_csv(OUTPUT_FILE, index=False)
    
#Merge genesetID to merged_fiveSpecies_component_summary
#
#Component summary rows:        56,370
#Component nodes with genesetid: 56,370
#check component nodes unique on componentID: True
#
#Merge results (component_summary left join genesetid):
#left_only         0
#right_only        0
#both          56370
#Name: merge_check, dtype: int64

if __name__ == '__main__':
    main()
