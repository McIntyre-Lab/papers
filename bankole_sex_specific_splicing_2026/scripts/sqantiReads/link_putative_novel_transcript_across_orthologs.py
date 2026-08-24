#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Check which orthologs across 5 species have novel isoforms
"""

import pandas as pd
from upsetplot import UpSet, from_indicators
import matplotlib.pyplot as plt

# === Load inputs ===
put_transcript_list = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list.csv'
ortholog_file = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/adp_ortholog_lists/ortholog_5species_list_02ndk.csv'

# Load data
transcript_df = pd.read_csv(put_transcript_list, index_col=0)
ortholog_df = pd.read_csv(ortholog_file)

# === Prepare subset of transcript_df with novel isoforms only ===
novel_pairs = transcript_df[['associated_gene', 'species', 'jxnHash']].drop_duplicates()

species_to_column = {
    'mel': 'dmel_geneID',
    'sim': 'dsim_geneID',
    'yak': 'dyak_geneID',
    'san': 'dsan_geneID',
    'ser': 'dser_geneID'
}

# === Create flags + collect jxnHash values for existing orthologs ===
for species, gene_col in species_to_column.items():
    subset = novel_pairs[novel_pairs['species'] == species]
    gene_to_jxnHash = (
        subset.groupby('associated_gene')['jxnHash']
        .apply(lambda hashes: ",".join(sorted(set(hashes))))
        .to_dict()
    )
    ortholog_df[f'flag_novel_transcript_in_{species}'] = ortholog_df[gene_col].astype(str).isin(gene_to_jxnHash).astype(int)
    ortholog_df[f'novel_jxnHash_{species}'] = ortholog_df[gene_col].astype(str).map(gene_to_jxnHash).fillna('')

# === Add novel genes not already in ortholog list ===
existing_gene_ids = set()
for species, gene_col in species_to_column.items():
    existing_gene_ids.update(ortholog_df[gene_col].dropna().astype(str).unique())

novel_genes_not_in_ortholog = novel_pairs[~novel_pairs['associated_gene'].astype(str).isin(existing_gene_ids)].copy()

# Create one-row-per-gene DataFrame with correct gene in correct species column
rows = []
for _, row in novel_genes_not_in_ortholog.iterrows():
    species = row['species']
    gene_id = row['associated_gene']
    jxnHash = row['jxnHash']
    
    new_row = {
        f'{species_to_column[species]}': gene_id,
        f'flag_novel_transcript_in_{species}': 1,
        f'novel_jxnHash_{species}': jxnHash,
        'flag_geneID_not_in_ortholog_list': 1
    }
    rows.append(new_row)

novel_gene_df = pd.DataFrame(rows)

# Add missing columns to match ortholog_df
for col in ortholog_df.columns:
    if col not in novel_gene_df.columns:
        novel_gene_df[col] = None

# Reorder columns, preserving custom flags
expected_cols = ortholog_df.columns.tolist()
if 'flag_geneID_not_in_ortholog_list' not in expected_cols:
    expected_cols.append('flag_geneID_not_in_ortholog_list')

novel_gene_df = novel_gene_df[[col for col in expected_cols if col in novel_gene_df.columns]]

# Combine new rows with ortholog_df
ortholog_augmented = pd.concat([ortholog_df, novel_gene_df], ignore_index=True)

# Fill NA in flags and jxnHash columns
for species in species_to_column:
    ortholog_augmented[f'flag_novel_transcript_in_{species}'] = ortholog_augmented[f'flag_novel_transcript_in_{species}'].fillna(0).astype(int)
    ortholog_augmented[f'novel_jxnHash_{species}'] = ortholog_augmented[f'novel_jxnHash_{species}'].fillna('')

# Fill or create flag_geneID_not_in_ortholog_list
if 'flag_geneID_not_in_ortholog_list' not in ortholog_augmented.columns:
    ortholog_augmented['flag_geneID_not_in_ortholog_list'] = 0
else:
    ortholog_augmented['flag_geneID_not_in_ortholog_list'] = ortholog_augmented['flag_geneID_not_in_ortholog_list'].fillna(0).astype(int)

# Summary flags
ortholog_augmented['flag_any_species_novel'] = ortholog_augmented[
    [f'flag_novel_transcript_in_{s}' for s in species_to_column]
].any(axis=1).astype(int)

ortholog_augmented['flag_novel_in_gt1_species'] = ortholog_augmented[
    [f'flag_novel_transcript_in_{s}' for s in species_to_column]
].sum(axis=1).ge(2).astype(int)

# Filter to only orthologs with at least one novel isoform
ortholog_with_novel = ortholog_augmented[ortholog_augmented['flag_any_species_novel'] == 1].copy()

#Merge in the ERP for putative novel transcript
species_to_tag = {
    'mel': 'dmel6',
    'sim': 'dsim2',
    'yak': 'dyak2',
    'san': 'dsan1',
    'ser': 'dser1'
}

for species, tag in species_to_tag.items():
    jxn_col = f'novel_jxnHash_{species}'
    erp_col = f'ERP_{species}'
    
    # Path to species-specific ERP file
    erp_file = f'/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/datafiles/datafile_jxnHash_{tag}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv'
    
    # Load only necessary columns
    try:
        erp_df = pd.read_csv(erp_file, usecols=[tag + '_jxnHash', 'ERP'])
    except FileNotFoundError:
        print(f"⚠️ File not found: {erp_file}")
        ortholog_with_novel[erp_col] = None
        continue

    # Rename column for merge clarity
    erp_df = erp_df.rename(columns={f'{tag}_jxnHash': 'jxnHash'})

    # Create mapping: jxnHash → ERP
    jxn_to_erp = erp_df.set_index('jxnHash')['ERP'].to_dict()
    
    # Map ERP values for all jxnHashes
    ortholog_with_novel[erp_col] = ortholog_with_novel[jxn_col].apply(
        lambda x: ",".join(
            [jxn_to_erp.get(j, '') for j in x.split(',') if j in jxn_to_erp]
        ) if pd.notna(x) and x != '' else None
    )

# === Output (commented out for testing) ===
ortholog_with_novel.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcripts_link_orthologs.csv', sep=',', index=False)

print("Done! Final dataset contains:")
print(f"- {ortholog_with_novel.shape[0]} ortholog groups with novel isoforms")
print(f"- {ortholog_augmented['flag_geneID_not_in_ortholog_list'].sum()} gene(s) added that were not in original ortholog table")
print(f"- {ortholog_augmented['flag_novel_in_gt1_species'].sum()} gene(s) with a putative novel trascnript in at least 2 species")

#Cretae upset plot
# Extract only the novel transcript flag columns
novel_flag_cols = [
    'flag_novel_transcript_in_mel',
    'flag_novel_transcript_in_sim',
    'flag_novel_transcript_in_yak',
    'flag_novel_transcript_in_san',
    'flag_novel_transcript_in_ser'
]

# Convert to boolean (if not already)
ortholog_with_novel[novel_flag_cols] = ortholog_with_novel[novel_flag_cols].astype(bool)

# Create a multi-index indicator DataFrame for UpSet plot
upset_data = from_indicators(novel_flag_cols, ortholog_with_novel)

# Generate the UpSet plot
plt.figure(figsize=(10, 6))
UpSet(upset_data, show_counts=True).plot()
plt.tight_layout()
plt.savefig("/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/upset_plot_novel_transcripts_link_orthologs.png")
