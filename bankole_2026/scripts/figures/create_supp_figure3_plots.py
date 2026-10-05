#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 12:46:50 2026

@author: mgaran
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PROJ    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{PROJ}/Figures"
SUMMARY_DIR = f"{PROJ}/zenodo/summary_files"
ROZ = "/TB14/TB14/mgaran"

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
    df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_{species}.csv",
        low_memory=False
    )

    multi_exon_genes = set(
        df.loc[
            (df['num_ujc'] - df['num_MonoExon']) >= 1,
            'geneID'
        ]
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
plt.savefig(f"{FIGURES_DIR}/supp_figure3_panelB.svg", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
print(f"Panel B saved to: {FIGURES_DIR}/supp_figure3_panelB.svg")



# PANEL C: transcript annotation proportions by observed non-ISM transcript count

category_order = ['1', '2', '3', '4', '5', '6+']
transcript_types = ['anno_UJC', 'anno_ERP', 'novel_ERP']

panel_c_data = {}

for species in species_list:
    df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_{species}.csv",
        low_memory=False
    )

    # Number of observed UJCs that are not ISM transcript fragments
    df['num_nonISM_UJC'] = (
        df['num_anno_ujc']
        + df['num_ujc_w_anno_erp']
        + df['num_ujc_w_novel_erp']
    )

    # Bin genes by number of observed non-ISM transcripts
    df['num_transcript_bin'] = df['num_nonISM_UJC'].apply(
        lambda x: '6+' if x >= 6 else str(int(x)) if x >= 1 else None
    )

    panel_c_data[species] = {}

    for category in category_order:
        df_cat = df[df['num_transcript_bin'] == category]

        # Sum transcript counts across genes in this transcript-count bin
        n_anno_UJC = df_cat['num_anno_ujc'].sum()
        n_anno_ERP = df_cat['num_ujc_w_anno_erp'].sum()
        n_novel_ERP = df_cat['num_ujc_w_novel_erp'].sum()

        # Denominator includes only observed non-ISM UJCs
        total_nonISM = n_anno_UJC + n_anno_ERP + n_novel_ERP

        panel_c_data[species][category] = {
            'anno_UJC': n_anno_UJC / total_nonISM if total_nonISM > 0 else 0,
            'anno_ERP': n_anno_ERP / total_nonISM if total_nonISM > 0 else 0,
            'novel_ERP': n_novel_ERP / total_nonISM if total_nonISM > 0 else 0,
            'total_nonISM': total_nonISM
        }

# Plot Panel C
fig, ax = plt.subplots(figsize=(16, 8))

x = np.arange(len(category_order))
bar_width = 0.12
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for i, species in enumerate(species_list):
    x_pos = x + offset_positions[i] * bar_width

    anno_ujc = np.array([
        panel_c_data[species][cat]['anno_UJC']
        for cat in category_order
    ])
    anno_erp = np.array([
        panel_c_data[species][cat]['anno_ERP']
        for cat in category_order
    ])
    novel_erp = np.array([
        panel_c_data[species][cat]['novel_ERP']
        for cat in category_order
    ])

    ax.bar(
        x_pos, anno_ujc, bar_width,
        color=dark[species],
        edgecolor='black',
        linewidth=0.5
    )
    ax.bar(
        x_pos, anno_erp, bar_width,
        bottom=anno_ujc,
        color=medium[species],
        edgecolor='black',
        linewidth=0.5
    )
    ax.bar(
        x_pos, novel_erp, bar_width,
        bottom=anno_ujc + anno_erp,
        color=light[species],
        edgecolor='black',
        linewidth=0.5
    )

ax.set_xlabel('Number of observed non-ISM transcripts per gene', fontsize=12, fontweight='bold')
ax.set_ylabel('Proportion of transcripts', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(category_order, fontsize=11)
ax.tick_params(axis='y', labelsize=11)
ax.set_ylim(0, 1.0)
ax.grid(axis='y', alpha=0.3, linestyle='--')

ax.legend(
    handles=[
        Patch(facecolor=dark['dmel'], edgecolor='black', label='Annotated transcript'),
        Patch(facecolor=medium['dmel'], edgecolor='black', label='Annotated exon pattern'),
        Patch(facecolor=light['dmel'], edgecolor='black', label='Novel exon pattern')
    ],
    loc='upper right',
    fontsize=9
)

plt.tight_layout()
plt.savefig(
    f"{FIGURES_DIR}/supp_figure3_panelC.svg",
    dpi=300,
    bbox_inches='tight'
)
plt.show()
plt.close()

print(f"Panel C saved to: {FIGURES_DIR}/supp_figure3_panelC.svg")


# PANEL D: AS comparison for genes with >=2 min10read non-ISM transcripts

panel_d_specs = [
    (['alt_donor_acceptor', 'alt_IR'],
     ['Alt. donor/acceptor', 'Alt. intron retention'],
     'Donor/acceptor and intron retention'),
    (['alt_5p_ER', 'alt_3p_ER', 'alt_ERSkip'],
     ['Alt. start (5′)', 'Alt. end (3′)', 'Alt. skip'],
     'Alternative start, end, and skip'),
]

plot_data = {}
for species in species_list:
    gene_df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_min10reads_{species}.csv",
        low_memory=False
    )
    as_df = pd.read_csv(
        f"{ROZ}/gene_summary_from_AS_analysis_min10reads_{species}.csv",
        low_memory=False
    ).set_index('geneID')

    multi_transcript = (gene_df['num_ujc'] - gene_df['num_ism_ujc']) >= 2
    multi_exon = gene_df['numExon_GM'] >= 2
    genes_plotted_df = gene_df[multi_transcript & multi_exon].copy()
    total = len(genes_plotted_df)
    print(f"Panel D - {species}: {total:,} genes in denominator")

    plot_data[species] = {}
    for alt_col in ['alt_donor_acceptor', 'alt_IR', 'alt_5p_ER', 'alt_3p_ER', 'alt_ERSkip']:
        tier_values = genes_plotted_df['geneID'].map(as_df[alt_col])
        plot_data[species][alt_col] = {
            'UJCanno': (tier_values == 'anno_UJC').sum() / total,
            'ERPanno': (tier_values == 'anno_ERP').sum() / total,
            'ERPnovel': (tier_values == 'novel_ERP').sum() / total,
        }

fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)
bar_width = 0.16

for ax, (alt_cols, alt_labels, subtitle) in zip(axes, panel_d_specs):
    x = np.arange(len(alt_cols))

    for i, species in enumerate(species_list):
        x_pos = x + i * bar_width
        ujc_vals = np.array([plot_data[species][col]['UJCanno'] for col in alt_cols])
        erp_a_vals = np.array([plot_data[species][col]['ERPanno'] for col in alt_cols])
        erp_n_vals = np.array([plot_data[species][col]['ERPnovel'] for col in alt_cols])

        ax.bar(x_pos, ujc_vals, bar_width, color=dark[species], edgecolor='black', linewidth=0.5)
        ax.bar(x_pos, erp_a_vals, bar_width, bottom=ujc_vals, color=medium[species], edgecolor='black', linewidth=0.5)
        ax.bar(x_pos, erp_n_vals, bar_width, bottom=ujc_vals + erp_a_vals, color=light[species], edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Alternative splicing category', fontsize=12)
    ax.set_ylabel('Proportion of genes', fontsize=10)
    ax.set_title(subtitle, fontsize=14)
    ax.set_xticks(x + (len(species_list) - 1) * bar_width / 2)
    ax.set_xticklabels(alt_labels, rotation=0)
    ax.set_ylim(0, 1.05)
    ax.legend(
        handles=[
            Patch(facecolor=dark['dmel'], edgecolor='black', label='UJCanno'),
            Patch(facecolor=medium['dmel'], edgecolor='black', label='ERPanno'),
            Patch(facecolor=light['dmel'], edgecolor='black', label='ERPnovel'),
        ],
        loc='upper right', fontsize=9
    )

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/as_comparison_min10reads.svg", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
print(f"Panel D saved to: {FIGURES_DIR}/as_comparison_min10reads.svg")
