#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 13:35:10 2026

@author: mgaran
"""

import pandas as pd
import numpy as np

BASE      = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
OUTDIR    = f"{BASE}/zenodo/summary_files"
ANNOT_DIR = f"{BASE}/zenodo"
DATA_DIR  = f"{BASE}/zenodo/datafiles"
NODE_PATH = f"{ANNOT_DIR}/fiveSpecies_supporting_files/network_files/component_map_by_node.csv"

SPECIES_GENOME_MAP = {
    'dmel': 'dmel6', 'dsim': 'dsim2', 'dser': 'dser1',
    'dsan': 'dsan1', 'dyak': 'dyak2'
}

# define the 2 groups for effect size and sexBias columns
G1 = 'F'
G2 = 'M'
CNT_PREFIX = 'rawCnts'

def calculate_component_bias_effectsize(df, g1, g2):
    """Calculate component_bias and component_effectSize from data aggregates.

    component_bias follows the column-definition logic exactly:
      'B' if > 1 female-biased jxnHash AND > 1 male-biased jxnHash
      Else g1 if > 1 female-biased jxnHash AND 0 male-biased jxnHash
      Else g2 if > 1 male-biased jxnHash AND 0 female-biased jxnHash
      Else 'U' if the number of analyzable jxnHash is > 1
      Else 'N'
    """
    df = df.copy()
    df['component_effectSize'] = np.where(
        (df['num_tested'] > 0) & (df['num_positive_effectSize'] == df['num_tested']), g1,
        np.where(
            (df['num_tested'] > 0) & (df['num_positive_effectSize'] == 0), g2,
            np.where((df['num_tested'] > 0), 'B', 'N')
        )
    )
    df['component_bias'] = np.where(
        (df[f'num_UJC_bias_{g1}'] > 1) & (df[f'num_UJC_bias_{g2}'] > 1), 'B',
        np.where(
            (df[f'num_UJC_bias_{g1}'] > 1) & (df[f'num_UJC_bias_{g2}'] == 0), g1,
            np.where(
                (df[f'num_UJC_bias_{g2}'] > 1) & (df[f'num_UJC_bias_{g1}'] == 0), g2,
                np.where((df['num_UJC_analyzable'] > 1), 'U', 'N')
            )
        )
    )
    return df

def process_species(species, genome, node_map, g1, g2):
    log_file = f"{OUTDIR}/log_files/create_component_summaries_{species}.log"
    print(f"Processing {species}...")

    # Load annotation file to get componentIDs for each jxnHash
    anno_df = pd.read_csv(
        f"{ANNOT_DIR}/fiveSpecies_{genome}_full_annotation.csv",
        low_memory=False
    )
    # rename columns for consistency
    anno_df = anno_df.rename(columns={
        'component_id': 'componentID',
        f'{genome}_jxnHash': 'jxnHash'
    })

    # flag annotated UJC with 1 exon
    anno_df['flag_monoExon_ujc'] = (anno_df['numExon_ERP'] == 1).astype(int)

    # aggregate jxnHash level columns in annotation file by componentID
    anno_comp = anno_df.groupby('componentID').agg(
        num_UJC=('jxnHash', 'nunique'),
        num_ERPp=('ERP_plus', 'nunique'),
        num_geneID=('geneID', 'nunique'),
        flag_monoExon=('flag_monoExon_ujc', 'min'),
    ).reset_index()

    # simple component: exactly 1 geneID and 1 ERP_plus (matches column definition)
    anno_comp['flag_simple'] = (
        (anno_comp['num_geneID'] == 1) & (anno_comp['num_ERPp'] == 1)
    ).astype(int)

    # merge geneID_concat onto annotation component summary
    geneID_concat = (
        anno_df.groupby('componentID')['geneID']
        .apply(lambda x: '|'.join(sorted(set(x.dropna()))))
        .reset_index(name='geneID_concat')
    )

    anno_comp = pd.merge(
        anno_comp, geneID_concat,
        on='componentID', how='outer', indicator='merge_check'
    )
    
    with open(log_file, 'w') as log:
        log.write(f"check that geneID_concat is unique on componentID before merging column to component_summary: {geneID_concat['componentID'].is_unique}")
        log.write(f"Merge geneID_concat to annotation component summary for {species}:\n")
        log.write(f"{anno_comp['merge_check'].value_counts(dropna=False).sort_index()}\n\n")
    anno_comp = anno_comp[anno_comp['merge_check'] == 'both'].drop('merge_check', axis=1)

    # Load datafile_jxnHash
    jxnHash_data = pd.read_csv(
        f"{DATA_DIR}/datafile_jxnHash_{genome}.csv",
        low_memory=False
    )
    # rename genome-specific jxnHash column to generic jxnHash
    jxnHash_data = jxnHash_data.rename(columns={f'{genome}_jxnHash': 'jxnHash'})
    if 'effect_size_Equal' not in jxnHash_data.columns:
        jxnHash_data['effect_size_Equal'] = np.nan
    # sum raw read count columns to get group and total read sums
    total_cnt_cols = [col for col in jxnHash_data.columns if col.startswith(f'{CNT_PREFIX}_')]
    jxnHash_data['sumReads'] = jxnHash_data[total_cnt_cols].sum(axis=1)
    for group in [g1, g2]:
        group_cols = [col for col in total_cnt_cols if f'_{group}_' in col]
        jxnHash_data[f'sumReads_{group}'] = jxnHash_data[group_cols].sum(axis=1) if group_cols else 0

    # merge jxnHash_data to node map on jxnHash
    node_map_sp = node_map[node_map['source'] == genome][
        ['jxnHash', 'geneID', 'componentID']
    ].copy()

    # drop geneID from the datafile (geneID comes from node_map_sp)
    jxnHash_data = jxnHash_data.drop(columns=['geneID'])

    jxnHash_data_w_compID = pd.merge(
        jxnHash_data,
        node_map_sp,
        on=['jxnHash'],
        how='outer',
        indicator='merge_check'
    )
    
    with open(log_file, 'a') as log:
        log.write(f"check datafile_jxnHash is unique on jxnHash: {jxnHash_data['jxnHash'].is_unique}\n")
        log.write(f"count unique jxnHash in datafile_jxnHash: {jxnHash_data['jxnHash'].nunique()}\n")
        log.write(f"check component map file is unique on jxnHash: {node_map_sp['jxnHash'].is_unique}\n")
        log.write(f"count unique jxnHash in component map file: {node_map_sp['jxnHash'].nunique()}\n")
        log.write(f"Merge node map to datafile_jxnHash for {species} (on jxnHash):\n")
        log.write(f"{jxnHash_data_w_compID['merge_check'].value_counts(dropna=False).sort_index()}\n\n")
        
    jxnHash_data_w_compID = jxnHash_data_w_compID[jxnHash_data_w_compID['merge_check'] == 'both'].drop('merge_check', axis=1)
    jxnHash_data_w_compID['componentID'] = jxnHash_data_w_compID['componentID'].astype(int)

    # per-junction flags
    jxnHash_data_w_compID[f'flag_bias_{g1}'] = (jxnHash_data_w_compID['sexClass'] == f'{g1}_bias').astype(int)
    jxnHash_data_w_compID[f'flag_bias_{g2}'] = (jxnHash_data_w_compID['sexClass'] == f'{g2}_bias').astype(int)
    jxnHash_data_w_compID['flag_tested']     = jxnHash_data_w_compID['sexClass'].isin(['unbiased', 'F_bias', 'M_bias']).astype(int)
    jxnHash_data_w_compID['flag_pos_es']     = (jxnHash_data_w_compID['effect_size_Equal'] > 0).astype(int)

    data_comp = jxnHash_data_w_compID.groupby('componentID').agg(
        sumReads=('sumReads', 'sum'),
        **{f'sumReads_{g1}': (f'sumReads_{g1}', 'sum')},
        **{f'sumReads_{g2}': (f'sumReads_{g2}', 'sum')},
        num_UJC_w_data=('jxnHash', 'nunique'),
        num_ERPp_w_data=('ERP_plus', 'nunique'),
        num_geneID_w_data=('geneID', 'nunique'),
        **{f'num_UJC_bias_{g1}': (f'flag_bias_{g1}', 'sum')},
        **{f'num_UJC_bias_{g2}': (f'flag_bias_{g2}', 'sum')},
        num_UJC_analyzable=('flag_analyzable', 'sum'),
        num_tested=('flag_tested', 'sum'),
        num_positive_effectSize=('flag_pos_es', 'sum'),
    ).reset_index()

    data_comp = calculate_component_bias_effectsize(data_comp, g1, g2)
    data_comp = data_comp.drop(['num_tested', 'num_positive_effectSize'], axis=1)

    # merge annotation-based and data-based component summaries
    # outer merge keeping all annotation components even if no data matched
    comp_df = pd.merge(
        anno_comp, data_comp,
        on='componentID', how='outer', indicator='merge_check'
    )
    with open(log_file, 'a') as log:
        log.write(f"Merge data-based component summary to annotation-based component summary for {species}:\n")
        log.write(f"{comp_df['merge_check'].value_counts(dropna=False).sort_index()}\n\n")
    comp_df = comp_df[comp_df['merge_check'].isin(['left_only', 'both'])].drop('merge_check', axis=1)

    with open(log_file, 'a') as log:
        log.write(f"Final component summary for {species}:\n")
        log.write(f"  Total components (annotation): {len(anno_comp):,}\n")
        log.write(f"  Components with data:          {len(data_comp):,}\n")
        log.write(f"  Components in merged output:   {len(comp_df):,}\n")

    # order columns
    cols = [
        'componentID',
        'flag_simple',
        'flag_monoExon',
        'geneID_concat',
        'num_UJC',
        'num_ERPp',
        'num_geneID',
        'num_UJC_w_data',
        'num_ERPp_w_data',
        'num_geneID_w_data',
        'sumReads',
        f'sumReads_{g1}',
        f'sumReads_{g2}',
        f'num_UJC_bias_{g1}',
        f'num_UJC_bias_{g2}',
        'num_UJC_analyzable',
        'component_bias',
    ]
    comp_df = comp_df[[c for c in cols if c in comp_df.columns]]

    output_file = f"{OUTDIR}/component_summary_{species}.csv"
    comp_df.to_csv(output_file, index=False)
    print(f"  Saved to: {output_file}")

def main():
    node_map = pd.read_csv(NODE_PATH, low_memory=False)
    node_map = node_map.rename(columns={'jxnhash': 'jxnHash', 'component_id': 'componentID'})

    for species, genome in SPECIES_GENOME_MAP.items():
        process_species(species, genome, node_map, G1, G2)


if __name__ == '__main__':
    main()
