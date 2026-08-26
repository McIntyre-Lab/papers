#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE_DIR    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures"
SUMMARY_DIR = f"{BASE_DIR}/zenodo/summary_files"

species_list = ['dmel', 'dsim', 'dsan', 'dyak', 'dser']

# Three-tier color palette (dark -> medium -> light)
dark = {'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A', 'dyak': '#717273', 'dser': '#825CA6'}
medium = {'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7', 'dyak': '#B8B9BA', 'dser': '#C3A9D3'}
light = {'dmel': '#F4D8A8', 'dsim': '#C8E2F8', 'dsan': '#C6E9E5', 'dyak': '#D9DADB', 'dser': '#E1D8F0'}

# PANEL B: multi-exon genes by non-fragment-transcript count and annotation tier

category_order = ['1', '2', '3', '4', '5', '6+']
tier_order     = ['anno_UJC', 'anno_ERP', 'novel_ERP']

panel_b_data = {}

for species in species_list:
    # Full-data summary defines which genes are multi-exon
    df_full = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_{species}.csv",
        low_memory=False
    )

    multi_exon_genes = set(
        df_full.loc[
            (df_full['num_ujc'] - df_full['num_MonoExon']) >= 1,
            'geneID'
        ]
    )

    # Min10 summary provides the transcript counts used for the bins
    df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_min10reads_{species}.csv",
        low_memory=False
    )

    # Keep only genes identified as multi-exon from the full-data summary
    df = df[df['geneID'].isin(multi_exon_genes)].copy()

    print(f"Panel B - {species}: {len(multi_exon_genes):,} multi-exon genes")
    # Non-fragment transcript count, binned (6+ for >= 6)
    df['num_nonISM_UJC'] = df['num_ujc'] - df['num_ism_ujc']
    df['num_transcript_bin'] = df['num_nonISM_UJC'].apply(
        lambda x: '6+' if x >= 6 else str(int(x)) if not pd.isna(x) else None
    )

    # Highest annotation tier present (hierarchical priority)
    df['hierarchy_tier'] = np.where(
        df['num_anno_ujc'] >= 1, 'anno_UJC',
        np.where(df['num_ujc_w_anno_erp'] >= 1, 'anno_ERP',
        np.where(df['num_ujc_w_novel_erp'] >= 1, 'novel_ERP', None))
    )

    panel_b_data[species] = {}
    for cat in category_order:
        df_cat = df[df['num_transcript_bin'] == cat]
        panel_b_data[species][cat] = {tier: (df_cat['hierarchy_tier'] == tier).sum() for tier in tier_order}

# Plot Panel B
fig, ax = plt.subplots(figsize=(16, 8))
x                = np.arange(len(category_order))
bar_width        = 0.12
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2
color_by_tier    = {'anno_UJC': dark, 'anno_ERP': medium, 'novel_ERP': light}

for i, species in enumerate(species_list):
    x_pos  = x + offset_positions[i] * bar_width
    bottom = np.zeros(len(category_order))
    for tier in tier_order:
        counts = [panel_b_data[species][cat][tier] for cat in category_order]
        ax.bar(x_pos, counts, bar_width, bottom=bottom,
               color=color_by_tier[tier][species], edgecolor='black', linewidth=0.5)
        bottom += np.array(counts)

ax.set_xlabel('num_nonISM_UJC per gene  (= num_ujc - num_ism_ujc)', fontsize=12, fontweight='bold')
ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
ax.set_title(
    'Number of multi-exon genes by annotation-status hierarchy, grouped by num_nonISM_UJC count\n'
    '(Dark = num_anno_ujc >= 1  |  Medium = num_ujc_w_anno_erp >= 1  |  Light = num_ujc_w_novel_erp >= 1)',
    fontsize=11, fontweight='bold'
)
ax.set_xticks(x)
ax.set_xticklabels(category_order, fontsize=11)
ax.tick_params(axis='y', labelsize=11)
ax.grid(axis='y', alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure3_panelB.svg", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
print(f"Panel B saved to: {FIGURES_DIR}/figure3_panelB.svg")

# PANELS C AND D: proportion of multiple-transcript (>=10 reads) genes with AS events

all_alt_columns = ['alt_donor_acceptor', 'alt_IR', 'alt_5p_ER', 'alt_3p_ER', 'alt_ERSkip']
plot_data = {col: {} for col in all_alt_columns}

for species in species_list:
    df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_min10reads_{species}.csv", low_memory=False)

    # AS event tiers from the min10 AS table (one row per geneID)
    as_df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_min10reads_from_AS_analysis_{species}.csv", low_memory=False
    ).set_index('geneID')
    for col in all_alt_columns:
        df[col] = df['geneID'].map(as_df[col])

    # Denominator: genes with multiple observed transcripts (>= 2 non-fragment)
    df['num_nonISM_UJC'] = df['num_ujc'] - df['num_ism_ujc']
    filtered_df = df[df['num_nonISM_UJC'] >= 2]
    total_count = len(filtered_df)
    print(f"Panels C/D - {species}: {total_count} multiple-transcript (>=10 reads) genes (denominator)")
    if total_count == 0:
        continue

    for col in all_alt_columns:
        n_ujc       = (filtered_df[col] == 'anno_UJC').sum()
        n_erp_anno  = (filtered_df[col] == 'anno_ERP').sum()
        n_erp_novel = (filtered_df[col] == 'novel_ERP').sum()
        plot_data[col][species] = {
            'anno_ujc':     n_ujc       / total_count,
            'anno_erp':     n_erp_anno  / total_count,
            'novel_erp':    n_erp_novel / total_count,
            'count_total':  n_ujc + n_erp_anno + n_erp_novel,
        }

legend_patches = [
    Patch(facecolor=dark['dmel'],   edgecolor='black', label='anno_ujc'),
    Patch(facecolor=medium['dmel'], edgecolor='black', label='anno_erp'),
    Patch(facecolor=light['dmel'],  edgecolor='black', label='novel_erp'),
]

# Panel C: alt_donor_acceptor and alt_IR
panel_c_cols   = ['alt_donor_acceptor', 'alt_IR']
panel_c_labels = ['Alt. donor/acceptor', 'Alt. intron retention']

fig, ax = plt.subplots(figsize=(9, 6))
x = np.arange(len(panel_c_cols)); bar_width = 0.16
for i, species in enumerate(species_list):
    x_pos      = x + i * bar_width
    ujc_vals   = np.array([plot_data[col][species]['anno_ujc']  for col in panel_c_cols])
    erp_a_vals = np.array([plot_data[col][species]['anno_erp']  for col in panel_c_cols])
    erp_n_vals = np.array([plot_data[col][species]['novel_erp'] for col in panel_c_cols])
    ax.bar(x_pos, ujc_vals,   bar_width, color=dark[species],   edgecolor='black', linewidth=0.5)
    ax.bar(x_pos, erp_a_vals, bar_width, bottom=ujc_vals, color=medium[species], edgecolor='black', linewidth=0.5)
    ax.bar(x_pos, erp_n_vals, bar_width, bottom=ujc_vals + erp_a_vals, color=light[species], edgecolor='black', linewidth=0.5)
ax.set_xlabel('Alternative splicing category', fontsize=12)
ax.set_ylabel('Proportion of genes', fontsize=10)
ax.set_title('Proportion of genes by annotation level of alternative splicing\n(donor/acceptor, intron retention)', fontsize=14)
ax.set_xticks(x + (len(species_list) - 1) * bar_width / 2)
ax.set_xticklabels(panel_c_labels, rotation=0)
ax.set_ylim(0, 1.05)
ax.legend(handles=legend_patches, loc='upper right', fontsize=9)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure3_panelC.svg", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
print(f"Panel C saved to: {FIGURES_DIR}/figure3_panelC.svg")

# Panel D: alt_5p_ER, alt_3p_ER, alt_ERSkip
panel_d_cols   = ['alt_5p_ER', 'alt_3p_ER', 'alt_ERSkip']
panel_d_labels = ['Alt. start (5p)', 'Alt. end (3p)', 'Alt. skip (ERSkip)']

fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(panel_d_cols)); bar_width = 0.16
for i, species in enumerate(species_list):
    x_pos      = x + i * bar_width
    ujc_vals   = np.array([plot_data[col][species]['anno_ujc']  for col in panel_d_cols])
    erp_a_vals = np.array([plot_data[col][species]['anno_erp']  for col in panel_d_cols])
    erp_n_vals = np.array([plot_data[col][species]['novel_erp'] for col in panel_d_cols])
    ax.bar(x_pos, ujc_vals,   bar_width, color=dark[species],   edgecolor='black', linewidth=0.5)
    ax.bar(x_pos, erp_a_vals, bar_width, bottom=ujc_vals, color=medium[species], edgecolor='black', linewidth=0.5)
    ax.bar(x_pos, erp_n_vals, bar_width, bottom=ujc_vals + erp_a_vals, color=light[species], edgecolor='black', linewidth=0.5)
ax.set_xlabel('Alternative splicing category', fontsize=12)
ax.set_ylabel('Proportion of genes', fontsize=10)
ax.set_title('Proportion of genes by annotation level of alternative splicing (5p, 3p, skip)', fontsize=14)
ax.set_xticks(x + (len(species_list) - 1) * bar_width / 2)
ax.set_xticklabels(panel_d_labels, rotation=0)
ax.set_ylim(0, 1.05)
ax.legend(handles=legend_patches, loc='upper right', fontsize=9)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure3_panelD.svg", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
print(f"Panel D saved to: {FIGURES_DIR}/figure3_panelD.svg")
