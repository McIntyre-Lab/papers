#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create all Figure 4 plots."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch

BASE_DIR    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures"
SUMMARY_DIR   = f"{BASE_DIR}/zenodo/summary_files"
TABLES_DIR    = f"{BASE_DIR}/Tables"
MERGED_ORTHOLOG_FILE = f"{TABLES_DIR}/gene_sexBias_merged_on_genesetid.csv"
DUP_DIR       = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/duplication_vs_splicing/analysis_adp"

species_list = ['dmel', 'dsim', 'dsan', 'dyak', 'dser']
sex_bias_species_list = ['dmel', 'dsim', 'dser']
analyzed_categories = ['F', 'M', 'B', 'unbiased']

dark   = {'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A', 'dyak': '#717273', 'dser': '#825CA6'}
medium = {'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7', 'dyak': '#B8B9BA', 'dser': '#C3A9D3'}
light  = {'dmel': '#F4D8A8', 'dsim': '#C8E2F8', 'dsan': '#C6E9E5', 'dyak': '#D9DADB', 'dser': '#E1D8F0'}
sexbias_color = {'F': '#FF0000', 'M': '#0000FF', 'B': '#A020F0'}

# PANEL A: proportion of evaluated genes with sex bias by transcript count
category_order = ['1', '2', '3', '4', '5', '6+']
panel_a_data = {}

for species in sex_bias_species_list:
    df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{species}.csv", low_memory=False)

    # Number of analyzable non-ISM UJCs per gene
    df['num_nonISM_UJC_analyzable'] = (
        df['num_anno_ujc_analyzable']
        + df['num_ujc_w_anno_erp_analyzable']
        + df['num_ujc_w_novel_erp_analyzable']
    )

    # Bin genes by number of analyzable non-ISM UJCs
    df['numTranscript_bin'] = df['num_nonISM_UJC_analyzable'].apply(
        lambda x: '6+' if x >= 6 else str(int(x)) if x >= 1 else None
    )

    panel_a_data[species] = {}

    for category in category_order:
        df_category = df[df['numTranscript_bin'] == category]

        # Denominator: all genes evaluated for sex bias in this transcript-count bin
        evaluated = df_category['sexBias_anno_status'] != 'not_evaluated'
        denominator = evaluated.sum()

        n_anno_UJC = (
            evaluated
            & (df_category['sexBias_anno_status'] == 'anno_UJC_bias')
        ).sum()

        n_anno_ERP = (
            evaluated
            & (df_category['sexBias_anno_status'] == 'anno_ERP_bias')
        ).sum()

        n_novel_ERP = (
            evaluated
            & (df_category['sexBias_anno_status'] == 'novel_ERP_bias')
        ).sum()

        n_ISM = (
            evaluated
            & (df_category['sexBias_anno_status'] == 'ISM_bias')
        ).sum()

        panel_a_data[species][category] = {
            'anno_UJC': n_anno_UJC / denominator if denominator > 0 else 0,
            'anno_ERP': n_anno_ERP / denominator if denominator > 0 else 0,
            'novel_ERP': n_novel_ERP / denominator if denominator > 0 else 0,
            'ISM': n_ISM / denominator if denominator > 0 else 0,
            'count': denominator
        }

    n_evaluated = (df['sexBias_anno_status'] != 'not_evaluated').sum()
    print(f"{species}: {n_evaluated} genes evaluated for sex bias")

fig, ax = plt.subplots(figsize=(14, 7))

x = np.arange(len(category_order))
bar_width = 0.15
offset_positions = (
    np.arange(len(sex_bias_species_list))
    - (len(sex_bias_species_list) - 1) / 2
)

for i, species in enumerate(sex_bias_species_list):
    xp = x + offset_positions[i] * bar_width

    anno_ujc = np.array([
        panel_a_data[species][cat]['anno_UJC']
        for cat in category_order
    ])
    anno_erp = np.array([
        panel_a_data[species][cat]['anno_ERP']
        for cat in category_order
    ])
    novel_erp = np.array([
        panel_a_data[species][cat]['novel_ERP']
        for cat in category_order
    ])
    ism = np.array([
        panel_a_data[species][cat]['ISM']
        for cat in category_order
    ])

    ax.bar(
        xp, anno_ujc, bar_width,
        color=dark[species], edgecolor='black', linewidth=0.5
    )
    ax.bar(
        xp, anno_erp, bar_width,
        bottom=anno_ujc,
        color=medium[species], edgecolor='black', linewidth=0.5
    )
    ax.bar(
        xp, novel_erp, bar_width,
        bottom=anno_ujc + anno_erp,
        color=light[species], edgecolor='black', linewidth=0.5
    )
    ax.bar(
        xp, ism, bar_width,
        bottom=anno_ujc + anno_erp + novel_erp,
        color='white', edgecolor='black', linewidth=0.8
    )

    total_height = anno_ujc + anno_erp + novel_erp + ism
    counts = [
        panel_a_data[species][cat]['count']
        for cat in category_order
    ]

    for xpos, height, count in zip(xp, total_height, counts):
        ax.text(
            xpos, height + 0.005, str(count),
            ha='center', va='bottom', fontsize=6, rotation=90
        )

ax.set_xlabel('Number of analyzable non-ISM transcripts per gene', fontsize=12)
ax.set_ylabel('Proportion of genes', fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(category_order)
ax.set_ylim(0, 0.9)
ax.grid(axis='y', alpha=0.3, linestyle='--')

ax.legend(
    handles=[
        Patch(
            facecolor=dark['dmel'],
            edgecolor='black',
            label='Annotated transcript'
        ),
        Patch(
            facecolor=medium['dmel'],
            edgecolor='black',
            label='Annotated exon pattern'
        ),
        Patch(
            facecolor=light['dmel'],
            edgecolor='black',
            label='Novel exon pattern'
        ),
        Patch(
            facecolor='white',
            edgecolor='black',
            label='Transcript fragment'
        )
    ],
    loc='upper left',
    fontsize=9
)

plt.tight_layout()
plt.savefig(
    f"{FIGURES_DIR}/figure4_panelA.svg",
    dpi=300,
    bbox_inches='tight'
)
plt.close()


# PANEL B: Gene copy number vs non-ISM transcript count (excluding fragments)
# Top row: Genes with sex bias
# Bottom row: Genes without statistical evidence for sex bias
rng = np.random.default_rng(seed=42)
panel_b = [
    ('dmel', f"{DUP_DIR}/dmel_singleton_vs_dup.csv", {'dmel': 'geneID'}, {}, 'geneID', 'num_FBGNs', 0.15, 'Melanogaster'),
    ('dsim', f"{DUP_DIR}/dsim_5species_Fbgn_2_Hahn_GLEANR.csv", {}, {'geneID': 'dsim_FBgn'}, 'dsim_FBgn', 'num_GLEANRIDs_sim', 0.001, 'Simulans'),
]

# Row specs: Top = with sex bias, Bottom = without sex bias (per legend)
row_specs = [
    ('w', 'Genes with sex bias'),
    ('wo', 'Genes without statistical evidence for sex bias')
]

fig, axes = plt.subplots(2, 4, figsize=(18, 12), sharey=True, gridspec_kw={'width_ratios': [2, 1, 2, 1], 'wspace': 0.08})

for col, (sp, dup_path, dup_rename, data_rename, merge_on, xcol, jit, name) in enumerate(panel_b):
    dup = pd.read_csv(dup_path).rename(columns=dup_rename)
    data = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{sp}.csv").rename(columns=data_rename)
    # Exclude transcript fragments (ISM fragments excluded by summing non-ISM ERPp categories)
    data['num_nonISM_ERPp'] = data['num_ERPp_w_anno_ujc'] + data['num_ERPp_w_anno_erp'] + data['num_ERPp_w_novel_erp']
    merged = pd.merge(dup, data, on=merge_on, how='inner')
    x = merged[xcol].to_numpy(dtype=float)
    y = merged['num_nonISM_ERPp'].to_numpy(dtype=float)
    mask = merged['sexBias'].isin(['M', 'F', 'B']).to_numpy()
    
    left_idx = col * 2
    right_idx = col * 2 + 1
    
    for r, (sub, tag) in enumerate(row_specs):
        m = mask if sub == 'w' else ~mask
        ax_l = axes[r, left_idx]   # Left subplot: 0 to 50
        ax_r = axes[r, right_idx]  # Right subplot: 150 to 250
        
        x_val = x[m]
        y_val = y[m]
        
        m_l = x_val <= 50
        m_r = x_val >= 150
        
        ax_l.scatter(x_val[m_l] + rng.uniform(-jit, jit, m_l.sum()), y_val[m_l], s=10, alpha=0.4, c=dark[sp], edgecolors='none', linewidths=0)
        ax_r.scatter(x_val[m_r] + rng.uniform(-jit, jit, m_r.sum()), y_val[m_r], s=10, alpha=0.4, c=dark[sp], edgecolors='none', linewidths=0)
        
        ax_l.set_xlim(0, 50)
        ax_r.set_xlim(150, 250)
        ax_l.set_ylim(0, 350)
        ax_r.set_ylim(0, 350)
        
        ax_l.spines['right'].set_visible(False)
        ax_r.spines['left'].set_visible(False)
        ax_l.yaxis.tick_left()
        ax_r.yaxis.tick_right()
        ax_r.tick_params(labelleft=False)
        
        if r == 1:
            ax_l.set_xlabel("Number of genes in gene family")
            ax_r.set_xlabel("")
        if left_idx == 0 and r == 0:
            ax_l.set_ylabel("number of ERP_plus without ISM reads")
            
        ax_l.set_title(f"{name} — {tag}", loc='left', fontsize=11, fontweight='bold')
        
        # Add diagonal break slash marks (//)
        d = 0.015
        kwargs = dict(transform=ax_l.transAxes, color='k', clip_on=False, linewidth=1)
        ax_l.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax_l.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
        
        kwargs.update(transform=ax_r.transAxes)
        ax_r.plot((-d, +d), (-d, +d), **kwargs)
        ax_r.plot((-d, +d), (1 - d, 1 + d), **kwargs)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure4_panelB.svg", dpi=300, bbox_inches='tight')
plt.close()


# PANEL C: proportion of analyzed genes with sex bias (per species + conserved)
sexbias_plot_data = {}
for species in sex_bias_species_list:
    df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{species}.csv", low_memory=False)
    analyzable = df['sexBias'].isin(analyzed_categories)
    n = analyzable.sum()
    sexbias_plot_data[species] = {
        'M_prop': (analyzable & (df['sexBias'] == 'M')).sum() / n,
        'F_prop': (analyzable & (df['sexBias'] == 'F')).sum() / n,
        'B_prop': (analyzable & (df['sexBias'] == 'B') & ((df['num_ujc'] - df['num_ism_ujc']) >= 2)).sum() / n,
    }

merged_df = pd.read_csv(MERGED_ORTHOLOG_FILE, low_memory=False)
for sp in sex_bias_species_list:
    merged_df[f'{sp}_multi'] = (merged_df[f'{sp}_num_ujc'] - merged_df[f'{sp}_num_ism_ujc']) >= 2
comparisons = [
    ('dmel and dser',        ['dmel_sexBias', 'dser_sexBias'],             ['dmel_multi', 'dser_multi'],                  'dmel_geneID'),
    ('dmel and dsim',        ['dmel_sexBias', 'dsim_sexBias'],             ['dmel_multi', 'dsim_multi'],                  'dmel_geneID'),
    ('dser and dsim',        ['dser_sexBias', 'dsim_sexBias'],             ['dser_multi', 'dsim_multi'],                  'dsim_geneID'),
    ('dmel, dsim, and dser', ['dmel_sexBias', 'dsim_sexBias', 'dser_sexBias'], ['dmel_multi', 'dsim_multi', 'dser_multi'], 'dmel_geneID'),
]

gap = 0.3
bars = [(i + 1, sexbias_plot_data[sp]) for i, sp in enumerate(sex_bias_species_list)]
comp_x0 = len(sex_bias_species_list) + 1 + gap
comp_labels = []
for j, (label, cols, mcols, gid) in enumerate(comparisons):
    denom = merged_df[cols].isin(analyzed_categories).all(axis=1)
    multi = merged_df[mcols].all(axis=1)
    n = merged_df.loc[denom, gid].nunique()
    props = {
        'M_prop': merged_df.loc[denom & (merged_df[cols] == 'M').all(axis=1), gid].nunique() / n,
        'F_prop': merged_df.loc[denom & (merged_df[cols] == 'F').all(axis=1), gid].nunique() / n,
        'B_prop': merged_df.loc[denom & (merged_df[cols] == 'B').all(axis=1) & multi, gid].nunique() / n,
    }
    bars.append((comp_x0 + j, props))
    comp_labels.append(label)

fig, ax = plt.subplots(figsize=(16, 7))
bar_width = 0.5
for xp, props in bars:
    ax.bar(xp, props['M_prop'], bar_width, color=sexbias_color['M'], edgecolor='black', linewidth=0.5)
    ax.bar(xp, props['F_prop'], bar_width, bottom=props['M_prop'], color=sexbias_color['F'], edgecolor='black', linewidth=0.5)
    ax.bar(xp, props['B_prop'], bar_width, bottom=props['M_prop'] + props['F_prop'], color=sexbias_color['B'], edgecolor='black', linewidth=0.5)
ax.axvline(x=(len(sex_bias_species_list) + comp_x0) / 2, color='gray', linestyle='--', linewidth=1, alpha=0.6)
ax.set_xticks([xp for xp, _ in bars])
ax.set_xticklabels(sex_bias_species_list + comp_labels, fontsize=9)
ax.set_xlabel('Species', fontsize=12); ax.set_ylabel('Proportion of Analyzable Genes', fontsize=12)
ax.set_title('Gene sex bias\n(Purple = both-sex biased, Red = Female biased, Blue = Male biased)\n'
             'Left: per-species proportions | Right: conserved sex bias across species pairs/trio', fontsize=12)
ax.set_ylim(0, 1.0); ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.legend(handles=[mpatches.Patch(facecolor=sexbias_color['B'], edgecolor='black', label='Both (B)'),
                   mpatches.Patch(facecolor=sexbias_color['F'], edgecolor='black', label='Female (F)'),
                   mpatches.Patch(facecolor=sexbias_color['M'], edgecolor='black', label='Male (M)')],
          loc='upper right', fontsize=10)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure4_panelC.svg", dpi=150, bbox_inches='tight')
plt.close()

# PANELS E: AS between M-biased and F-biased UJCs among B-biased genes
mvf_plots = [
    (['alt_donor_acceptor_FvM', 'alt_IR_FvM'], ['Alt. donor/acceptor FvM', 'Alt. intron retention FvM'],
     'donor/acceptor, intron retention', 'figure4_panelE_donor_acceptor_IR.svg', (9, 6)),
    (['alt_5p_ER_FvM', 'alt_3p_ER_FvM', 'alt_ERSkip_FvM'], ['Alt. start (5\u2032) FvM', 'Alt. end (3\u2032) FvM', 'Alt. skip FvM'],
     '5\u2032, 3\u2032, skip', 'figure4_panelE_5p_3p_skip.svg', (10, 6)),
]
for alt_cols, alt_labels, subtitle, outname, figsize in mvf_plots:
    plot_data = {lab: {} for lab in alt_labels}
    for species in sex_bias_species_list:
        df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{species}.csv", low_memory=False)
        as_df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_from_AS_analysis_{species}.csv", low_memory=False).set_index('geneID')
        for c in alt_cols:
            df[c] = df['geneID'].map(as_df[c])
        fdf = df[((df['num_ujc'] - df['num_ism_ujc']) >= 2) & (df['sexBias'] == 'B')]
        tot = len(fdf)
        for c, lab in zip(alt_cols, alt_labels):
            plot_data[lab][species] = {
                'UJCanno': (fdf[c] == 'anno_UJC').sum() / tot,
                'ERPanno': (fdf[c] == 'anno_ERP').sum() / tot,
                'ERPnovel': (fdf[c] == 'novel_ERP').sum() / tot,
            }
    fig, ax = plt.subplots(figsize=figsize)
    bar_width = 0.16
    xs = np.arange(len(alt_labels))
    for i, species in enumerate(sex_bias_species_list):
        xp = xs + i * bar_width
        u = np.array([plot_data[l][species]['UJCanno'] for l in alt_labels])
        a = np.array([plot_data[l][species]['ERPanno'] for l in alt_labels])
        nov = np.array([plot_data[l][species]['ERPnovel'] for l in alt_labels])
        ax.bar(xp, u, bar_width, color=dark[species], edgecolor='black', linewidth=0.5)
        ax.bar(xp, a, bar_width, bottom=u, color=medium[species], edgecolor='black', linewidth=0.5)
        ax.bar(xp, nov, bar_width, bottom=u + a, color=light[species], edgecolor='black', linewidth=0.5)
    ax.set_xlabel('Alternative splicing category', fontsize=12); ax.set_ylabel('Proportion of genes', fontsize=12)
    ax.set_title(f'Proportion of genes with AS between M-biased and F-biased UJCs\n({subtitle})', fontsize=14)
    ax.set_xticks(xs + (len(sex_bias_species_list) - 1) * bar_width / 2)
    ax.set_xticklabels(alt_labels, rotation=0); ax.set_ylim(0, 1.05)
    ax.legend(handles=[Patch(facecolor=dark['dmel'], edgecolor='black', label='UJCanno'),
                       Patch(facecolor=medium['dmel'], edgecolor='black', label='ERPanno'),
                       Patch(facecolor=light['dmel'], edgecolor='black', label='ERPnovel')], loc='upper right', fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/{outname}", dpi=300, bbox_inches='tight')
    plt.close()