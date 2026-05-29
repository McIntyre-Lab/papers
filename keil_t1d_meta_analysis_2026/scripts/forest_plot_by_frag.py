#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 00:14:27 2025

@author: nkeil
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import matplotlib.patches as patches


path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/data_4_forest_plots'
graphs_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/graphs/forest_plots_by_frag'

os.makedirs(graphs_path, exist_ok=True)

# ----------------------------
# Collect study files
# ----------------------------
all_files = os.listdir(path)
study_files = [f for f in all_files if f.endswith('_study.csv')]

# ----------------------------
# Helpers
# ----------------------------
def sanitize_filename(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

def safe_fragment_label(feature_id: str) -> str:
    # remove gene prefix (ENSG...:)
    if ':' in feature_id:
        return ':'.join(feature_id.split(':')[1:])
    return feature_id

# ----------------------------
# Plot a single fragment (F vs M)
# ----------------------------
def plot_one_fragment(df_female_row, df_male_row, gene_name, cell_type, frag_label, outdir):
    """
    df_female_row and df_male_row are either a pandas Series (one row) or None if missing for that sex.
    """
    # Determine y-limits from whichever rows are present
    ci_lowers, ci_uppers = [], []
    if df_female_row is not None:
        ci_lowers.append(df_female_row['ci_lower'])
        ci_uppers.append(df_female_row['ci_upper'])
    if df_male_row is not None:
        ci_lowers.append(df_male_row['ci_lower'])
        ci_uppers.append(df_male_row['ci_upper'])

    if not ci_lowers or not ci_uppers:
        return

    ymin = min(ci_lowers) - 0.25
    ymax = max(ci_uppers) + 0.25

    # Figure
    plt.figure(figsize=(6, 6))
    ax = plt.gca()

    # x positions for F and M
    x_f, x_m = 0, 1
    fixed_width = 0.25
    base_height = 1.0

    # Legend
    plt.plot([], [], 's', color='red', alpha=1.0, label='Female (exon)')
    plt.plot([], [], 's', color='red', alpha=0.5, linestyle='dashed', label='Female (intron)')
    plt.plot([], [], 'o', color='blue', alpha=1.0, label='Male (exon)')
    plt.plot([], [], 'o', color='blue', alpha=0.5, linestyle='dashed', label='Male (intron)')

    # Female marker
    if df_female_row is not None:
        alpha_f = 0.5 if int(df_female_row['ef_ir_flag']) == 1 else 1.0
        linestyle_f = 'dashed' if int(df_female_row['ef_ir_flag']) == 1 else 'solid'
        height_f = base_height * (float(df_female_row['weight']) / 100.0)
        rect = patches.Rectangle(
            (x_f - fixed_width / 2, float(df_female_row['effect_size']) - height_f / 2),
            fixed_width, height_f, linewidth=1,
            edgecolor='red', facecolor='red', alpha=alpha_f, linestyle=linestyle_f
        )
        ax.add_patch(rect)
        plt.errorbar(
            x_f, float(df_female_row['effect_size']),
            yerr=[
                [float(df_female_row['effect_size']) - float(df_female_row['ci_lower'])],
                [float(df_female_row['ci_upper']) - float(df_female_row['effect_size'])]
            ],
            fmt='none', ecolor='red', alpha=alpha_f, capsize=5
        )

    # Male marker
    if df_male_row is not None:
        alpha_m = 0.5 if int(df_male_row['ef_ir_flag']) == 1 else 1.0
        linestyle_m = 'dashed' if int(df_male_row['ef_ir_flag']) == 1 else 'solid'
        height_m = base_height * (float(df_male_row['weight']) / 100.0)
        oval = patches.Ellipse(
            (x_m, float(df_male_row['effect_size'])),
            fixed_width, height_m, linewidth=1,
            edgecolor='blue', facecolor='blue', alpha=alpha_m, linestyle=linestyle_m
        )
        ax.add_patch(oval)
        plt.errorbar(
            x_m, float(df_male_row['effect_size']),
            yerr=[
                [float(df_male_row['effect_size']) - float(df_male_row['ci_lower'])],
                [float(df_male_row['ci_upper']) - float(df_male_row['effect_size'])]
            ],
            fmt='none', ecolor='blue', alpha=alpha_m, capsize=5
        )

    # Styling
    plt.axhline(0, color='grey', linestyle='--', alpha=0.5)
    plt.xlim(-0.5, 1.5)
    plt.ylim(ymin, ymax)
    plt.xticks([x_f, x_m], ['Female', 'Male'])
    plt.ylabel('Effect Size', labelpad=12)
    plt.title(f'{gene_name} ({cell_type}) — {frag_label}')
    plt.grid(True, alpha=0.4)
    plt.legend(loc='upper right')
    plt.tight_layout()

    # Save
    fname = f'ef_sz_{sanitize_filename(gene_name)}_{sanitize_filename(cell_type)}_{sanitize_filename(frag_label)}.png'
    plt.savefig(os.path.join(outdir, fname), dpi=300)
    plt.close()

# ----------------------------
# Main: loop files, then fragments
# ----------------------------
for file in study_files:
    gene_name_match = re.search(r'ENSG\d+', file)
    if not gene_name_match:
        continue
    gene_name = gene_name_match.group(0)
    cell_type_match = re.search(r'CD[48]', file)
    cell_type = cell_type_match.group(0) if cell_type_match else 'Unknown'

    df = pd.read_csv(os.path.join(path, file))
    df['frag_label'] = df['featureID'].apply(safe_fragment_label)

    if 'ef_start' in df.columns:
        df = df.sort_values(by=['ef_start'], ascending=True)

    df_female = df[df['sex'] == 'F']
    df_male = df[df['sex'] == 'M']

    unique_feats = pd.unique(df['featureID'])

    for feat in unique_feats:
        frag_label = df.loc[df['featureID'] == feat, 'frag_label'].iloc[0]

        f_row = df_female[df_female['featureID'] == feat]
        m_row = df_male[df_male['featureID'] == feat]

        f_row = f_row.iloc[0] if len(f_row) > 0 else None
        m_row = m_row.iloc[0] if len(m_row) > 0 else None

        if f_row is None and m_row is None:
            continue

        plot_one_fragment(f_row, m_row, gene_name, cell_type, frag_label, graphs_path)