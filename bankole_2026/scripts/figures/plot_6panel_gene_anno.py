#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:38:24 2026

@author: mgaran
"""

import os
import gc
import colorsys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

from plot_6panel_utils import (
    ALL_SPECIES_DEFAULT,
    SPECIES_NAMES,
    load_gtf_coords_by_jxnhash,
    flip_coords_to_display,
    load_components_from_list,
    create_component_alignment,
    sort_components_by_erp,
    plot_component_panel,
    compute_xlim,
    apply_x_formatting,
    compute_shared_xlim,
    parse_species_paths,
    build_index_csv,
)

# ---------------------------------------------------------------------------
# Species color scheme (gene_anno only)
# ---------------------------------------------------------------------------

SPECIES_COLORS = {
    'dmel6': '#966729',
    'dsim2': '#3F78C1',
    'dsan1': '#28827A',
    'dyak2': '#717273',
    'dser1': '#825CA6',
}

# Flag columns indicating a UJC was present in the original annotation.
# For dsim2, either flag == 1 counts.
ORIGINAL_ANNO_FLAGS = {
    'dmel6': ['flag_dmel650_2_dmel6_ujc'],
    'dsim2': ['flag_dsim202_2_dsim2_ujc', 'flag_dsimWXD_2_dsim2_ujc'],
    'dyak2': ['flag_dyak21_2_dyak2_ujc'],
    'dsan1': ['flag_dsan11_2_dsan1_ujc'],
    'dser1': ['flag_dser11_2_dser1_ujc'],
}

JXNHASH_COL = {
    'dmel6': 'dmel6_jxnHash',
    'dsim2': 'dsim2_jxnHash',
    'dyak2': 'dyak2_jxnHash',
    'dsan1': 'dsan1_jxnHash',
    'dser1': 'dser1_jxnHash',
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def darken_hex(hex_color, factor=0.55):
    """Return a darkened version of a hex color."""
    hex_color = hex_color.lstrip('#')
    r, g, b = (int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0, l * factor)
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return '#{:02x}{:02x}{:02x}'.format(int(r2*255), int(g2*255), int(b2*255))


def load_original_anno_flags(full_anno_file, sp):
    """Return a set of jxnHashes that were in the original annotation for sp."""
    jxn_col   = JXNHASH_COL[sp]
    flag_cols = ORIGINAL_ANNO_FLAGS[sp]

    needed = [jxn_col] + flag_cols
    df = pd.read_csv(full_anno_file, low_memory=False,
                     usecols=lambda c: c in needed)

    if jxn_col not in df.columns:
        return set()

    present_flags = [f for f in flag_cols if f in df.columns]
    if not present_flags:
        return set()

    mask = df[present_flags].eq(1).any(axis=1)
    return set(df.loc[mask, jxn_col].dropna().unique())


def plot_anno_species_panel(ax, coords, positions, sp, max_pos,
                             original_jxns, show_xaxis=True, forced_xlim=None):
    """Plot annotated UJCs for one species using the species color scheme."""
    normal_color = SPECIES_COLORS[sp]
    dark_color   = darken_hex(normal_color)

    ax.set_title(SPECIES_NAMES[sp], fontweight='bold', fontsize=10)
    ax.set_yticks([])

    if not coords:
        ax.set_xlim(*(forced_xlim if forced_xlim else (0, 1)))
        ax.set_ylim(0.5, max_pos + 0.5)
        apply_x_formatting(ax, forced_xlim or (0, 1), 0, show_xaxis)
        ax.text(0.5, 0.5, 'no annotated UJCs', ha='center', va='center',
                fontsize=8, color='#999999', transform=ax.transAxes)
        return

    min_c, max_c = compute_xlim(coords, positions)
    pad  = 0.05 * (max_c - min_c) if max_c > min_c else 1000
    xlim = forced_xlim if forced_xlim else (min_c, max_c)

    for jxn, exons in coords.items():
        if jxn not in positions:
            continue
        y     = positions[jxn]
        color = dark_color if jxn in original_jxns else normal_color
        prev_end = None
        for start, end in exons:
            ax.add_patch(Rectangle(
                (start, y - 0.3), end - start, 0.6,
                facecolor=color, edgecolor='none', alpha=0.85,
            ))
            if prev_end is not None:
                ax.plot([prev_end, start], [y, y],
                        color=color, alpha=0.85, linewidth=1,
                        solid_capstyle='butt')
            prev_end = end

    apply_x_formatting(ax, xlim, pad if not forced_xlim else 0, show_xaxis)
    ax.set_ylim(0.5, max_pos + 0.5)


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def plot_gene_anno(component_file, component_list_file,
                   anno_gtf_files, full_anno_files,
                   gene_name, output_dir,
                   output_format='png'):
    """Create the annotation-only species-color 6-panel plot."""

    all_species  = ALL_SPECIES_DEFAULT
    comp_list_df = pd.read_csv(component_list_file)

    if 'component_id' not in comp_list_df.columns:
        raise ValueError("component_list file must contain 'component_id' column")

    if 'geneID' in comp_list_df.columns and 'source_species' in comp_list_df.columns:
        panel_title = (f"{comp_list_df['geneID'].iloc[0]}"
                       f":{comp_list_df['source_species'].iloc[0]}")
    else:
        panel_title = 'Component ID'

    component_ids = comp_list_df['component_id'].tolist()
    components, _ = load_components_from_list(component_file, component_ids, all_species)

    # ------------------------------------------------------------------
    # Load GTF coordinates per species
    # ------------------------------------------------------------------
    all_coords           = {}
    all_dominant_strands = {}

    for sp in all_species:
        sp_jxns = [j for cd in components.values() for j in cd[sp]]

        if not sp_jxns or sp not in anno_gtf_files:
            all_coords[sp]           = {}
            all_dominant_strands[sp] = None
            continue

        raw_coords, strands, _ = load_gtf_coords_by_jxnhash(
            anno_gtf_files[sp], sp_jxns, return_chr=True)

        strand_vals     = list(strands.values())
        dominant_strand = (max(set(strand_vals), key=strand_vals.count)
                           if strand_vals else '+')
        all_dominant_strands[sp] = dominant_strand
        all_coords[sp]           = flip_coords_to_display(raw_coords, dominant_strand)

    component_order = sort_components_by_erp(comp_list_df)

    master_positions, annotated_max_pos, component_positions = (
        create_component_alignment(components, all_species, all_coords,
                                   component_order=component_order))

    annotated_max_pos = max(annotated_max_pos, 5)

    # ------------------------------------------------------------------
    # Load original annotation flags per species
    # ------------------------------------------------------------------
    original_jxns = {}
    for sp in all_species:
        if sp in full_anno_files and os.path.exists(full_anno_files[sp]):
            original_jxns[sp] = load_original_anno_flags(full_anno_files[sp], sp)
        else:
            print(f'warning: full annotation file not found for {sp}')
            original_jxns[sp] = set()

    # ------------------------------------------------------------------
    # Compute shared x-limits per species
    # ------------------------------------------------------------------
    shared_xlim = {sp: compute_shared_xlim([all_coords[sp]]) for sp in all_species}

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    widths     = [0.4, 1, 1, 1, 1, 1]
    fig_height = max(annotated_max_pos * 0.3, 4)
    fig = plt.figure(figsize=(24, fig_height))
    gs  = fig.add_gridspec(1, 6, width_ratios=widths,
                           left=0.08, right=0.95, top=0.85, bottom=0.15,
                           wspace=0.25)
    axes = [fig.add_subplot(gs[0, i]) for i in range(6)]

    plot_component_panel(axes[0], component_positions, annotated_max_pos,
                         title=panel_title)

    for idx, sp in enumerate(all_species, start=1):
        plot_anno_species_panel(
            axes[idx],
            all_coords[sp],
            master_positions[sp],
            sp,
            annotated_max_pos,
            original_jxns[sp],
            show_xaxis=True,
            forced_xlim=shared_xlim.get(sp),
        )

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------
    legend_elements = []
    for sp in all_species:
        normal = SPECIES_COLORS[sp]
        dark   = darken_hex(normal)
        legend_elements.append(
            Patch(facecolor=normal, label=f'{SPECIES_NAMES[sp]} (cross-species)'))
        legend_elements.append(
            Patch(facecolor=dark,   label=f'{SPECIES_NAMES[sp]} (original annotation)'))

    fig.legend(handles=legend_elements, loc='lower center',
               ncol=5, fontsize=7, bbox_to_anchor=(0.5, -0.02))

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    os.makedirs(output_dir, exist_ok=True)
    out_file = f'{output_dir}/{gene_name}_anno_ujc_species_color.{output_format}'
    save_kw  = {} if output_format == 'svg' else {'dpi': 150}
    fig.savefig(out_file, bbox_inches='tight', **save_kw)
    print(f'saved: {out_file}')
    plt.close(fig)
    gc.collect()

    # ------------------------------------------------------------------
    # Save CSV lookup table
    # ------------------------------------------------------------------
    csv_rows = []
    for sp in all_species:
        if not all_coords[sp]:
            continue
        jxn_list_ordered = sorted(
            master_positions[sp].keys(),
            key=lambda j: master_positions[sp][j]
        )
        empty_reads = {j: 0 for j in jxn_list_ordered}
        sp_csv = build_index_csv(
            'UJCanno', sp, jxn_list_ordered,
            empty_reads, empty_reads, empty_reads
        )
        csv_rows.append(sp_csv)
    
    if csv_rows:
        csv_df = pd.concat(csv_rows, ignore_index=True)
        csv_file = out_file.replace(f'.{output_format}', '_index_lookup.csv')
        csv_df.to_csv(csv_file, index=False)
        print(f'saved: {csv_file}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Annotation-only species-color 6-panel plot — gene mode'
    )
    parser.add_argument('--component_file',  required=True)
    parser.add_argument('--component_list',  required=True)
    parser.add_argument('--gene_name',       required=True)
    parser.add_argument('--output_dir',      required=True)
    parser.add_argument('--anno_gtf_files',  nargs='+', required=True,
                        metavar='SPECIES:PATH')
    parser.add_argument('--full_anno_files', nargs='+', required=True,
                        metavar='SPECIES:PATH')
    parser.add_argument('--output_format',   choices=['svg', 'png', 'pdf'], default='png')
    args = parser.parse_args()

    plot_gene_anno(
        component_file      = args.component_file,
        component_list_file = args.component_list,
        anno_gtf_files      = parse_species_paths(args.anno_gtf_files),
        full_anno_files     = parse_species_paths(args.full_anno_files),
        gene_name           = args.gene_name,
        output_dir          = args.output_dir,
        output_format       = args.output_format,
    )