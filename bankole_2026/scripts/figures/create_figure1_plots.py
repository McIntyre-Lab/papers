#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import defaultdict

PROJ    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO  = f"{PROJ}/zenodo"
DATA    = f"{ZENODO}/datafiles"
ANNO    = f"{ZENODO}/fiveSpecies_supporting_files"
NETWORK = f"{ANNO}/network_files"
SUMMARY = f"{ZENODO}/summary_files"
OUT     = f"{PROJ}/Figures"

species_list = ['dmel', 'dsim', 'dyak', 'dsan', 'dser']

species_2_genome = {
    'dmel': 'dmel6', 'dsim': 'dsim2', 'dsan': 'dsan1',
    'dyak': 'dyak2', 'dser': 'dser1'
}

species_colors_detected = {
    'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A',
    'dyak': '#717273', 'dser': '#825CA6'
}
species_colors_undetected = {
    'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7',
    'dyak': '#B8B9BA', 'dser': '#C3A9D3'
}

###############################################################################
# PANEL C: Sankey diagram (component origin -> component annotation)
df = pd.read_csv(f"{PROJ}/Tables/cnt_origin_anno_xpress_byspecies.csv")
count_col = [c for c in df.columns if c.lower() in ('count', 'cnt', 'n')][0]
df = (
    df.groupby(['origin', 'comp_anno_in'], dropna=False)[count_col]
    .sum()
    .reset_index(name='COUNT')
)
 
# Build node labels with prefixes to keep them unique across sides
origins    = df['origin'].unique().tolist()
comp_annos = df['comp_anno_in'].unique().tolist()
all_labels = [f"Origin: {x}" for x in origins] + [f"Comp Anno: {x}" for x in comp_annos]
label_to_idx = {label: idx for idx, label in enumerate(all_labels)}
 
# Aggregate COUNT per (origin, comp_anno_in) pair
link_counts = defaultdict(int)
for _, row in df.iterrows():
    src = label_to_idx[f"Origin: {row['origin']}"]
    tgt = label_to_idx[f"Comp Anno: {row['comp_anno_in']}"]
    link_counts[(src, tgt)] += row['COUNT']
 
sources = [k[0] for k in link_counts]
targets = [k[1] for k in link_counts]
values  = list(link_counts.values())
 
fig = go.Figure(data=[go.Sankey(
    node=dict(pad=15, thickness=20, line=dict(color='black', width=0.5),
              label=all_labels, color='blue'),
    link=dict(source=sources, target=targets, value=values)
)])
fig.update_layout(title="Sankey Plot: Origin -> Comp Anno", font=dict(size=12),
                  height=800, width=1200)
fig.write_image(f"{OUT}/figure1_panelC_sankey_origin_2_anno.svg", format='svg', scale=4)
print(f"Sankey saved to: {OUT}/figure1_panelC_sankey_origin_2_anno.svg")

###############################################################################
# PANEL D top: genes by number of exons in gene model
species_results = {}

for species in species_list:

    print(f"Processing {species}...")

    anno_df = pd.read_csv(
        f"{ZENODO}/fiveSpecies_{species_2_genome[species]}_full_annotation.csv",
        low_memory=False
    )

    summary_df = pd.read_csv(
        f"{SUMMARY}/gene_summary_{species}.csv",
        usecols=['geneID'],
        low_memory=False
    )

    # One row per gene using the first ERP in the fiveSpecies annotation
    gene_df = anno_df.groupby('geneID')['ERP'].first().reset_index()

    # Flag detection status
    gene_df['flag_detected'] = gene_df['geneID'].isin(summary_df['geneID']).astype(int)

    # Calculate number of exons in the gene model
    gene_df['numExon_GM'] = gene_df['ERP'].str.len() - 2

    # Flag exon-count categories
    gene_df['flag_numExon_1'] = (gene_df['numExon_GM'] == 1).astype(int)
    gene_df['flag_numExon_2'] = (gene_df['numExon_GM'] == 2).astype(int)
    gene_df['flag_numExon_3'] = (gene_df['numExon_GM'] == 3).astype(int)
    gene_df['flag_numExon_4'] = (gene_df['numExon_GM'] == 4).astype(int)
    gene_df['flag_numExon_5plus'] = (gene_df['numExon_GM'] >= 5).astype(int)

    # Flag GFF5 genes
    gene_df['flag_GFF5'] = gene_df['geneID'].str.startswith('GFF').astype(int)

    species_results[species] = gene_df

# Plot conditions
gene_type_dict = {
    'GFF5': {'flag': 1, 'title': 'GFF5 genes by numExon_GM'},
    'nonGFF5': {'flag': 0, 'title': 'non-GFF5 genes by numExon_GM'}
}

detection_dict = {
    'detected': 1,
    'undetected': 0
}

exon_dict = {
    '1': 'flag_numExon_1',
    '2': 'flag_numExon_2',
    '3': 'flag_numExon_3',
    '4': 'flag_numExon_4',
    '5+': 'flag_numExon_5plus'
}

fig, axes = plt.subplots(1, 2, figsize=(20, 7))
x = np.arange(len(exon_dict))
bar_width = 0.15
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for ax_i, (gene_type, gene_type_info) in enumerate(gene_type_dict.items()):

    ax = axes[ax_i]

    for species_i, species in enumerate(species_list):

        gene_df = species_results[species]
        count_data = {}

        for detection_name, detection_flag in detection_dict.items():
            count_data[detection_name] = {}

            for exon_name, exon_flag in exon_dict.items():
                condition = (
                    (gene_df['flag_GFF5'] == gene_type_info['flag'])
                    & (gene_df['flag_detected'] == detection_flag)
                    & (gene_df[exon_flag] == 1)
                )
                count_data[detection_name][exon_name] = int(condition.sum())

        x_pos = x + offset_positions[species_i] * bar_width
        detected_vals = list(count_data['detected'].values())
        undetected_vals = list(count_data['undetected'].values())

        label_detected = f'{species} (detected)' if gene_type == 'GFF5' else ''
        label_undetected = f'{species} (undetected)' if gene_type == 'GFF5' else ''

        ax.bar(
            x_pos, detected_vals, bar_width,
            color=species_colors_detected[species],
            edgecolor='black', linewidth=0.5,
            label=label_detected
        )
        ax.bar(
            x_pos, undetected_vals, bar_width,
            bottom=detected_vals,
            color=species_colors_undetected[species],
            edgecolor='black', linewidth=0.5,
            label=label_undetected
        )

    ax.set_xlabel('Number of Exons in gene model', fontsize=11)
    ax.set_ylabel('Number of Genes', fontsize=11)
    ax.set_title(gene_type_info['title'], fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(exon_dict.keys(), fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUT}/figure1_panelD_genes_by_numExon.svg", bbox_inches='tight')
plt.close()

print(f"Saved: {OUT}/figure1_panelD_genes_by_numExon.svg")

###############################################################################
# PANEL D bottom: genes by transcript category

species_results = {}

for species in species_list:

    print(f"Processing {species}...")

    anno_df = pd.read_csv(
        f"{ZENODO}/fiveSpecies_{species_2_genome[species]}_full_annotation.csv",
        low_memory=False
    )

    summary_df = pd.read_csv(
        f"{SUMMARY}/gene_summary_{species}.csv",
        usecols=['geneID'],
        low_memory=False
    )

    # One row per gene using the first ERP in the fiveSpecies annotation
    gene_df = anno_df.groupby('geneID')['ERP'].first().reset_index()

    # Flag detection status
    gene_df['flag_detected'] = gene_df['geneID'].isin(summary_df['geneID']).astype(int)

    # Calculate number of exons in the gene model
    gene_df['numExon_GM'] = gene_df['ERP'].str.len() - 2

    # Flag exon-count categories
    gene_df['flag_numExon_1'] = (gene_df['numExon_GM'] == 1).astype(int)
    gene_df['flag_numExon_2'] = (gene_df['numExon_GM'] == 2).astype(int)
    gene_df['flag_numExon_3'] = (gene_df['numExon_GM'] == 3).astype(int)
    gene_df['flag_numExon_4'] = (gene_df['numExon_GM'] == 4).astype(int)
    gene_df['flag_numExon_5plus'] = (gene_df['numExon_GM'] >= 5).astype(int)

    # Flag GFF5 genes
    gene_df['flag_GFF5'] = gene_df['geneID'].str.startswith('GFF').astype(int)

    species_results[species] = gene_df

###############################################################################

# Plot conditions
gene_type_dict = {
    'GFF5': {'flag': 1, 'title': 'GFF5 genes by numExon_GM'},
    'nonGFF5': {'flag': 0, 'title': 'non-GFF5 genes by numExon_GM'}
}

detection_dict = {
    'detected': 1,
    'undetected': 0
}

exon_dict = {
    '1': 'flag_numExon_1',
    '2': 'flag_numExon_2',
    '3': 'flag_numExon_3',
    '4': 'flag_numExon_4',
    '5+': 'flag_numExon_5plus'
}

fig, axes = plt.subplots(1, 2, figsize=(20, 7))
x = np.arange(len(exon_dict))
bar_width = 0.15
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for ax_i, (gene_type, gene_type_info) in enumerate(gene_type_dict.items()):

    ax = axes[ax_i]

    for species_i, species in enumerate(species_list):

        gene_df = species_results[species]
        count_data = {}

        for detection_name, detection_flag in detection_dict.items():
            count_data[detection_name] = {}

            for exon_name, exon_flag in exon_dict.items():
                condition = (
                    (gene_df['flag_GFF5'] == gene_type_info['flag'])
                    & (gene_df['flag_detected'] == detection_flag)
                    & (gene_df[exon_flag] == 1)
                )
                count_data[detection_name][exon_name] = int(condition.sum())

        x_pos = x + offset_positions[species_i] * bar_width
        detected_vals = list(count_data['detected'].values())
        undetected_vals = list(count_data['undetected'].values())

        label_detected = f'{species} (detected)' if gene_type == 'GFF5' else ''
        label_undetected = f'{species} (undetected)' if gene_type == 'GFF5' else ''

        ax.bar(
            x_pos, detected_vals, bar_width,
            color=species_colors_detected[species],
            edgecolor='black', linewidth=0.5,
            label=label_detected
        )
        ax.bar(
            x_pos, undetected_vals, bar_width,
            bottom=detected_vals,
            color=species_colors_undetected[species],
            edgecolor='black', linewidth=0.5,
            label=label_undetected
        )

    ax.set_xlabel('Number of Exons in gene model', fontsize=11)
    ax.set_ylabel('Number of Genes', fontsize=11)
    ax.set_title(gene_type_info['title'], fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(exon_dict.keys(), fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUT}/figure1_panelD_genes_by_numExon.svg", bbox_inches='tight')
plt.close()

print(f"Saved: {OUT}/figure1_panelD_genes_by_numExon.svg")