#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 15:38:25 2026

@author: mgaran
"""

import os
import argparse
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import ScalarFormatter

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NO_DATA_COLOR = '#D3D3D3'

ALL_SPECIES_DEFAULT = ['dmel6', 'dsim2', 'dyak2', 'dsan1', 'dser1']

SPECIES_NAMES = {
    'dmel6': 'D. melanogaster',
    'dsim2': 'D. simulans',
    'dyak2': 'D. yakuba',
    'dsan1': 'D. santomea',
    'dser1': 'D. serrata',
}

DATA_NAMES = {
    'dmel6': 'dmel',
    'dsim2': 'dsim',
    'dyak2': 'dyak',
    'dsan1': 'dsan',
    'dser1': 'dser',
}


# ---------------------------------------------------------------------------
# CLI helper
# ---------------------------------------------------------------------------

def parse_species_paths(values):
    """Parse a list of 'species:path' strings into a dict."""
    result = {}
    for item in values:
        if ':' not in item:
            raise argparse.ArgumentTypeError(
                f'Entries must be "species:path", got: {item!r}')
        sp, path = item.split(':', 1)
        result[sp.strip()] = path.strip()
    return result


# ---------------------------------------------------------------------------
# X-limit computation
# ---------------------------------------------------------------------------

def compute_shared_xlim(coords_tiers):
    """Return a padded (min, max) xlim from the first non-empty coord dict.

    Parameters
    ----------
    coords_tiers : list of dict
        Priority-ordered list of jxnHash -> exon_list dicts.
        The first non-empty dict is used; the rest are fallbacks.

    Returns
    -------
    tuple or None
        (disp_min - pad, disp_max + pad), or None if all tiers are empty.
    """
    for coords in coords_tiers:
        exons = [e for exon_list in coords.values() for e in exon_list]
        if exons:
            disp_min = min(e[0] for e in exons)
            disp_max = max(e[1] for e in exons)
            pad = 0.05 * (disp_max - disp_min) if disp_max > disp_min else 1000
            return (disp_min - pad, disp_max + pad)
    return None


# ---------------------------------------------------------------------------
# GTF / coordinate helpers
# ---------------------------------------------------------------------------

def load_gtf_coords_by_jxnhash(gtf_file, jxn_hashes, return_chr=False):
    """Load exon coordinates (and optionally chromosome) from a GTF file."""
    df = pd.read_csv(
        gtf_file, sep='\t', comment='#', header=None,
        names=['chr', 'src', 'feat', 'start', 'end',
               'score', 'strand', 'frame', 'attr'],
        usecols=['chr', 'feat', 'start', 'end', 'strand', 'attr'],
        low_memory=False,
    )
    df = df[df['feat'] == 'exon']
    df['jxnHash'] = df['attr'].str.extract(r'transcript_id "([^"]+)"')
    df = df[df['jxnHash'].isin(jxn_hashes)]

    coords = (
        df.groupby('jxnHash')
        .apply(lambda x: sorted(zip(x['start'], x['end'])))
        .to_dict()
    )
    strands = df.groupby('jxnHash')['strand'].first().to_dict()

    if return_chr:
        chrs = df.groupby('jxnHash')['chr'].first().to_dict()
        return coords, strands, chrs

    return coords, strands


def auto_flip_by_strand(coords_dict, strands_dict):
    """Flip coordinates for minus-strand junctions so they display L→R."""
    flipped = {}
    for jxn, exons in coords_dict.items():
        if strands_dict.get(jxn, '+') == '-':
            flipped[jxn] = sorted(
                [(-end, -start) for start, end in exons], key=lambda x: x[0]
            )
        else:
            flipped[jxn] = exons
    return flipped


def flip_coords_to_display(raw_coords_dict, dominant_strand):
    """Flip all coords using the dominant strand direction."""
    if dominant_strand == '-':
        return {
            jxn: sorted([(-end, -start) for start, end in exons],
                        key=lambda x: x[0])
            for jxn, exons in raw_coords_dict.items()
        }
    return {jxn: list(exons) for jxn, exons in raw_coords_dict.items()}


def get_gene_region_from_raw_coords(raw_coords_dict):
    """Return the genomic (unflipped) min/max across all jxnHashes."""
    if not raw_coords_dict:
        return None, None
    all_exons = [e for exons in raw_coords_dict.values() for e in exons]
    return min(e[0] for e in all_exons), max(e[1] for e in all_exons)


# ---------------------------------------------------------------------------
# Component / datafile helpers
# ---------------------------------------------------------------------------

def load_components_from_list(component_file, component_list, all_species):
    """Load component → jxnHash mapping from the component map CSV."""
    df = pd.read_csv(component_file, low_memory=False)
    df = df.rename(columns={'jxnhash': 'jxnHash'})
    df = df[df['component_id'].isin(component_list)]

    components = {}
    comp_gene_ids = {sp: {} for sp in all_species}

    for comp_id, group in df.groupby('component_id'):
        components[comp_id] = {sp: [] for sp in all_species}
        for _, row in group.iterrows():
            if row['source'] in all_species:
                components[comp_id][row['source']].append(row['jxnHash'])
                if 'geneID' in row and pd.notna(row['geneID']):
                    comp_gene_ids[row['source']][row['jxnHash']] = row['geneID']

    return components, comp_gene_ids


def load_jxn_datafile(datafile_path, species, jxn_hashes):
    """Load only the rows for the given jxnHashes from a jxnHash datafile."""
    col_map = {
        'dmel6': ['rawCnts_dmel_F_rep4', 'rawCnts_dmel_F_rep5', 'rawCnts_dmel_F_rep6',
                  'rawCnts_dmel_M_rep4', 'rawCnts_dmel_M_rep5', 'rawCnts_dmel_M_rep6'],
        'dsim2': ['rawCnts_dsim_F_rep4', 'rawCnts_dsim_F_rep5', 'rawCnts_dsim_F_rep6',
                  'rawCnts_dsim_M_rep4', 'rawCnts_dsim_M_rep5', 'rawCnts_dsim_M_rep6'],
        'dsan1': ['rawCnts_dsan_F_rep1', 'rawCnts_dsan_M_rep1'],
        'dyak2': ['rawCnts_dyak_F_rep1', 'rawCnts_dyak_M_rep1'],
        'dser1': ['rawCnts_dser_F_rep4', 'rawCnts_dser_F_rep5', 'rawCnts_dser_F_rep6',
                  'rawCnts_dser_M_rep4', 'rawCnts_dser_M_rep5', 'rawCnts_dser_M_rep6'],
    }
    jxn_set = set(jxn_hashes)
    species_count_cols = set(col_map.get(species, []))
    
    # 1. Define the dynamic species hash name and add it to needed
    sp_hash = f"{species}_jxnHash"
    needed = {'jxnHash', sp_hash, 'effect_size_Equal', 'flag_ERPanno', 'geneID', 'ERP_plus'} | species_count_cols

    hits = []
    for chunk in pd.read_csv(datafile_path, low_memory=False, chunksize=50_000,
                              usecols=lambda c: c in needed):
        # 2. Swap out the hardcoded 'dmel6_jxnHash' with our dynamic variable
        if sp_hash in chunk.columns: chunk = chunk.rename(columns={sp_hash: 'jxnHash'})
        
        match = chunk[chunk['jxnHash'].isin(jxn_set)]
        if not match.empty:
            hits.append(match)
    return pd.concat(hits, ignore_index=True) if hits else pd.DataFrame(columns=list(needed))


def load_read_props_by_jxnhash(jxn_df, species, jxn_hashes):
    """Load read data and flags for a list of jxnHashes."""
    col_map = {
        'dmel6': ['rawCnts_dmel_F_rep4', 'rawCnts_dmel_F_rep5', 'rawCnts_dmel_F_rep6',
                  'rawCnts_dmel_M_rep4', 'rawCnts_dmel_M_rep5', 'rawCnts_dmel_M_rep6'],
        'dsim2': ['rawCnts_dsim_F_rep4', 'rawCnts_dsim_F_rep5', 'rawCnts_dsim_F_rep6',
                  'rawCnts_dsim_M_rep4', 'rawCnts_dsim_M_rep5', 'rawCnts_dsim_M_rep6'],
        'dsan1': ['rawCnts_dsan_F_rep1', 'rawCnts_dsan_M_rep1'],
        'dyak2': ['rawCnts_dyak_F_rep1', 'rawCnts_dyak_M_rep1'],
        'dser1': ['rawCnts_dser_F_rep4', 'rawCnts_dser_F_rep5', 'rawCnts_dser_F_rep6',
                  'rawCnts_dser_M_rep4', 'rawCnts_dser_M_rep5', 'rawCnts_dser_M_rep6'],
    }

    df = jxn_df[jxn_df['jxnHash'].isin(jxn_hashes)].copy()

    if df.empty:
        return {}, {}, {}, {}, {}, {}, {}

    cols = col_map.get(species, [])
    f_cols = [c for c in cols if '_F_' in c and c in df.columns]
    m_cols = [c for c in cols if '_M_' in c and c in df.columns]

    df['total_F'] = df[f_cols].sum(axis=1) if f_cols else 0
    df['total_M'] = df[m_cols].sum(axis=1) if m_cols else 0
    df['total'] = df['total_F'] + df['total_M']
    df['prop_F'] = (df['total_F'] / df['total'].replace(0, float('nan')))

    effect_sizes = {}
    if 'effect_size_Equal' in df.columns:
        effect_sizes = df.set_index('jxnHash')['effect_size_Equal'].to_dict()

    total_reads = df.set_index('jxnHash')['total'].to_dict()

    flag_erp_anno = {}
    if 'flag_ERPanno' in df.columns:
        flag_erp_anno = df.set_index('jxnHash')['flag_ERPanno'].to_dict()

    gene_ids = {}
    if 'geneID' in df.columns:
        gene_ids = df.set_index('jxnHash')['geneID'].to_dict()

    sumReads_f = df.set_index('jxnHash')['total_F'].to_dict()
    sumReads_m = df.set_index('jxnHash')['total_M'].to_dict()
    
    return (df.set_index('jxnHash')['prop_F'].to_dict(),
            effect_sizes,
            total_reads,
            flag_erp_anno,
            gene_ids,
            sumReads_f,
            sumReads_m)


def sort_components_by_erp(component_list_df):
    """Sort component IDs by ERP column (NaN last)."""
    if 'ERP' not in component_list_df.columns:
        raise ValueError('component_list file is missing the ERP column')

    df = component_list_df.copy()
    df['sort_key'] = df.apply(
        lambda r: (1, str(r['component_id'])) if pd.isna(r['ERP'])
        else (0, str(r['ERP'])),
        axis=1,
    )
    return df.sort_values('sort_key')['component_id'].tolist()


def create_component_alignment(components, all_species, all_coords,
                                component_order=None):
    """Assign shared y-positions across species based on component membership."""
    master_positions = {sp: {} for sp in all_species}
    component_positions = {}
    pos = 1

    if component_order is None:
        component_order = sorted(components.keys())

    for comp_id in component_order:
        comp_data = components[comp_id]
        filtered = {
            sp: [j for j in comp_data[sp] if j in all_coords.get(sp, {})]
            for sp in all_species
        }
        max_ujcs = max((len(v) for v in filtered.values()), default=0)
        if max_ujcs == 0:
            continue

        component_positions[comp_id] = pos
        for i in range(max_ujcs):
            for sp in all_species:
                if i < len(filtered[sp]):
                    master_positions[sp][filtered[sp][i]] = pos
            pos += 1

    return master_positions, pos - 1, component_positions


# ---------------------------------------------------------------------------
# Unannotated UJC discovery
# ---------------------------------------------------------------------------

def load_unannotated_jxns_in_region(data_gtf_file, gene_chr, gene_min, gene_max,
                                     known_jxns, gene_strand=None):
    """Find novel jxnHashes in a data GTF overlapping a gene's genomic region."""
    if not os.path.exists(data_gtf_file):
        return []

    df = pd.read_csv(
        data_gtf_file, sep='\t', comment='#', header=None,
        names=['chr', 'src', 'feat', 'start', 'end',
               'score', 'strand', 'frame', 'attr'],
        usecols=['chr', 'feat', 'start', 'end', 'strand', 'attr'],
        low_memory=False,
    )
    df = df[df['feat'] == 'exon']
    df = df[df['chr'] == gene_chr]
    df['jxnHash'] = df['attr'].str.extract(r'transcript_id "([^"]+)"')
    df = df[(df['start'] <= gene_max) & (df['end'] >= gene_min)]

    if gene_strand is not None:
        df = df[df['strand'] == gene_strand]

    return list(set(df['jxnHash'].dropna().unique()) - set(known_jxns))


# ---------------------------------------------------------------------------
# ERP effect-size helpers
# ---------------------------------------------------------------------------

def load_erp_effect_sizes_for_gene(erp_datafile, gene_id):
    """Load ERP_plus -> effect_size_Equal for all tested ERPs in a gene."""
    if not os.path.exists(erp_datafile):
        return {}

    df = pd.read_csv(erp_datafile, low_memory=False,
                     usecols=lambda c: c in {'geneID', 'ERP_plus', 'effect_size_Equal'})

    if 'geneID' not in df.columns:
        return {}

    df = df[df['geneID'] == gene_id]

    if 'effect_size_Equal' not in df.columns:
        return {}

    df = df[pd.notna(df['effect_size_Equal'])]

    if df.empty:
        return {}

    return df.set_index('ERP_plus')['effect_size_Equal'].to_dict()


def get_jxnhash_to_erp_map(jxn_df, gene_id):
    """Return jxnHash -> ERP_plus for all jxnHashes of a gene."""
    if 'geneID' not in jxn_df.columns or 'ERP_plus' not in jxn_df.columns:
        return {}
    df = jxn_df[jxn_df['geneID'] == gene_id]

    if df.empty:
        return {}

    return df.set_index('jxnHash')['ERP_plus'].to_dict()


def get_jxnhash_to_erp_map_for_jxns(jxn_df, jxn_hashes):
    """Return jxnHash -> ERP_plus for a specific list of jxnHashes."""
    if 'ERP_plus' not in jxn_df.columns or 'jxnHash' not in jxn_df.columns:
        return {}
    df = jxn_df[jxn_df['jxnHash'].isin(jxn_hashes)]

    if df.empty:
        return {}

    return df.set_index('jxnHash')['ERP_plus'].to_dict()


def build_jxn_erp_effect_sizes(jxn_datafile, erp_datafile, gene_id):
    """Build jxnHash -> ERP_plus effect_size for a gene."""
    jxn_to_erp = get_jxnhash_to_erp_map(jxn_datafile, gene_id)
    if not jxn_to_erp:
        return {}

    erp_to_es = load_erp_effect_sizes_for_gene(erp_datafile, gene_id)
    if not erp_to_es:
        return {}

    return {
        jxn: erp_to_es[erp]
        for jxn, erp in jxn_to_erp.items()
        if erp in erp_to_es
    }


# ---------------------------------------------------------------------------
# Core plotting helpers
# ---------------------------------------------------------------------------

def ujc_colors(prop_f, effect_size, colormap, normalize):
    """Return (fill_color, hatch_color, add_pattern, is_significant).

    Coloring scheme:
      fill        = effect_size_Equal  (blue→purple→red; gray if untested)
      hatch_color = prop_F colored when no effect size data
      pattern     = hatch added when no effect size data (gray fill)
    """
    prop_norm = Normalize(vmin=0, vmax=1)

    if pd.notna(effect_size):
        fc = colormap(normalize(max(-9, min(9, effect_size))))
    else:
        fc = NO_DATA_COLOR

    if pd.notna(prop_f):
        hatch_color = colormap(prop_norm(prop_f))
    else:
        hatch_color = '#888888'

    is_sig = (pd.notna(effect_size)
              and (effect_size > 2.78 or effect_size < -2.78))
    add_pat = (fc == NO_DATA_COLOR)

    return fc, hatch_color, add_pat, is_sig


def plot_ujc(ax, exons, y, fill_color, hatch_color='#888888', alpha=0.8,
             add_pattern=False, is_significant=False, panel_min=None,
             erp_effect_size=None, colormap=None, normalize=None):
    """Draw a single UJC (rectangles + intron lines) at y-position *y*.

    Markers (drawn at the left panel edge):
      - Yellow star  : UJC itself has a significant effect size (is_significant).
      - Filled circle: UJC has no significant effect size of its own, but belongs
                       to an ERP_plus that is significant (|erp_effect_size| > 2.78).
                       Colored by the ERP_plus effect size.
    Star takes priority — if is_significant, no circle is drawn.
    """
    prev_end = None
    for start, end in exons:
        ax.add_patch(Rectangle(
            (start, y - 0.3), end - start, 0.6,
            facecolor=fill_color,
            edgecolor=hatch_color if add_pattern else 'none',
            linewidth=0.3 if add_pattern else 0,
            hatch='///' if add_pattern else None,
            alpha=alpha,
        ))
        if prev_end is not None:
            ax.plot([prev_end, start], [y, y],
                    color=fill_color, alpha=alpha,
                    linewidth=1, solid_capstyle='butt')
        prev_end = end

    if panel_min is None:
        return

    marker_x = panel_min

    if is_significant:
        ax.plot(marker_x, y, marker='*', markersize=10, color='#FFFF00',
                markeredgecolor='black', markeredgewidth=0.5, zorder=10,
                clip_on=False)
    elif (erp_effect_size is not None
          and pd.notna(erp_effect_size)
          and abs(erp_effect_size) > 2.78
          and colormap is not None
          and normalize is not None):
        circle_color = colormap(normalize(max(-9, min(9, erp_effect_size))))
        ax.plot(marker_x, y, marker='o', markersize=8,
                color=circle_color,
                markeredgecolor='black', markeredgewidth=0.5, zorder=10,
                clip_on=False)


def compute_xlim(coords, positions):
    """Return (min_coord, max_coord) for UJCs that have a position."""
    plotted = [c for jxn, exons in coords.items()
               if jxn in positions for c in exons]
    if not plotted:
        plotted = [c for exons in coords.values() for c in exons]
    if not plotted:
        return 0, 1
    return min(c[0] for c in plotted), max(c[1] for c in plotted)


def apply_x_formatting(ax, xlim, pad, show_xaxis):
    ax.set_xlim(xlim[0] - pad, xlim[1] + pad)
    if show_xaxis:
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.ticklabel_format(style='plain', axis='x')
        ax.set_xlabel('Genomic Position', fontsize=9)
        ax.tick_params(axis='x', labelsize=8)
    else:
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', bottom=False)


def plot_annotated_panel(ax, coords, positions, props, effect_sizes,
                         species_name, max_pos,
                         colormap, normalize,
                         erp_es_by_jxn=None,
                         show_xaxis=False, forced_xlim=None):
    """Plot the annotated UJC row for one species."""
    ax.set_title(species_name, fontweight='bold', fontsize=10)
    ax.set_yticks([])

    if not coords:
        ax.set_xlim(*(forced_xlim if forced_xlim else (0, 1)))
        ax.set_ylim(0.5, max_pos + 0.5)
        apply_x_formatting(ax, forced_xlim or (0, 1), 0, show_xaxis)
        return

    min_c, max_c = compute_xlim(coords, positions)
    pad = 0.05 * (max_c - min_c) if max_c > min_c else 1000
    xlim = forced_xlim if forced_xlim else (min_c, max_c)
    left_edge = xlim[0] if forced_xlim else min_c - pad

    for jxn, exons in coords.items():
        if jxn not in positions:
            continue
        fc, hc, add_pat, is_sig = ujc_colors(
            props.get(jxn), effect_sizes.get(jxn), colormap, normalize)
        erp_es = (erp_es_by_jxn or {}).get(jxn)
        plot_ujc(ax, exons, positions[jxn], fc, hatch_color=hc,
                 add_pattern=add_pat, is_significant=is_sig,
                 panel_min=left_edge,
                 erp_effect_size=erp_es,
                 colormap=colormap, normalize=normalize)

    apply_x_formatting(ax, xlim, pad if not forced_xlim else 0, show_xaxis)
    ax.set_ylim(0.5, max_pos + 0.5)


def plot_novel_panel(ax, coords, props, effect_sizes,
                     species_name, colormap, normalize,
                     erp_es_by_jxn=None,
                     jxn_erp_group=None,
                     panel_label='novel',
                     show_xaxis=True, forced_xlim=None,
                     max_rows=None):
    """Plot a novel UJC row for one species.

    Parameters
    ----------
    max_rows : int or None
        If provided, ylim is set to this value so all species panels in the
        row have the same UJC height. Defaults to the local UJC count.

    Returns
    -------
    int  number of UJCs plotted
    """
    ax.set_title('')
    ax.set_yticks([])

    if not coords:
        ylim_max = max_rows if max_rows else 1
        ax.set_xlim(*(forced_xlim if forced_xlim else (0, 1)))
        ax.set_ylim(0.5, ylim_max + 0.5)
        ax.text(0.5, 0.5, f'no {panel_label} UJCs', ha='center', va='center',
                fontsize=8, color='#999999', transform=ax.transAxes)
        apply_x_formatting(ax, forced_xlim or (0, 1), 0, show_xaxis)
        return 0

    erp_map = jxn_erp_group or {}

    def erp_sort_key(jxn):
        erp = erp_map.get(jxn)
        if erp:
            return (0, str(erp))
        return (1, '')

    sorted_jxns = sorted(coords.keys(), key=erp_sort_key)
    n = len(sorted_jxns)
    positions = {jxn: i + 1 for i, jxn in enumerate(sorted_jxns)}

    min_c, max_c = compute_xlim(coords, positions)
    pad = 0.05 * (max_c - min_c) if max_c > min_c else 1000
    xlim = forced_xlim if forced_xlim else (min_c, max_c)
    left_edge = xlim[0] if forced_xlim else min_c - pad

    for jxn in sorted_jxns:
        exons = coords[jxn]
        y = positions[jxn]
        fc, hc, add_pat, is_sig = ujc_colors(
            props.get(jxn), effect_sizes.get(jxn), colormap, normalize)
        erp_es = (erp_es_by_jxn or {}).get(jxn)
        plot_ujc(ax, exons, y, fc, hatch_color=hc,
                 add_pattern=add_pat, is_significant=is_sig,
                 panel_min=left_edge,
                 erp_effect_size=erp_es,
                 colormap=colormap, normalize=normalize)

    apply_x_formatting(ax, xlim, pad if not forced_xlim else 0, show_xaxis)
    ylim_max = max_rows if max_rows else n
    ax.set_ylim(0.5, ylim_max + 0.5)

    return n


def plot_component_panel(ax, component_positions, max_pos, title='Component ID'):
    """Draw the leftmost component-ID label panel for the annotated row."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, max_pos + 0.5)
    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    for comp_id, pos in component_positions.items():
        ax.text(0.5, pos, str(comp_id), ha='center', va='center', fontsize=7)


def plot_novel_label_panel(ax, max_novel_rows, label='Novel UJCs'):
    """Draw the leftmost label panel for a novel UJC row with row indices."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0.5, max(max_novel_rows, 1) + 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for i in range(1, max_novel_rows + 1):
        ax.text(0.5, i, str(i), ha='center', va='center', fontsize=7)


def get_colormap():
    """Return (colormap, normalizer) for effect_size coloring."""
    colormap = LinearSegmentedColormap.from_list('mf', ['blue', 'purple', 'red'])
    normalize = Normalize(vmin=-9, vmax=9)
    return colormap, normalize

def build_index_csv(annotation_status, species, jxn_list, sumReads_f, sumReads_m, sumReads):
    """Build CSV table for one annotation_status x species panel.
    
    Parameters:
        annotation_status: 'UJCanno', 'ERPanno', or 'ERPnovel'
        species: species code (dmel6, dsim2, etc.)
        jxn_list: ordered list of jxnHashes as they appear in plot
        sumReads_f: dict of jxnHash -> female read count
        sumReads_m: dict of jxnHash -> male read count
        sumReads: dict of jxnHash -> total read count
    
    Returns: pd.DataFrame with 7 columns
    """
    rows = []
    for idx, jxn in enumerate(jxn_list, start=1):
        rows.append({
            'annotation_status': annotation_status,
            'species': species,
            'index': idx,
            'jxnHash': jxn,
            'sumReads_f': sumReads_f.get(jxn, 0),
            'sumReads_m': sumReads_m.get(jxn, 0),
            'sumReads': sumReads.get(jxn, 0),
        })
    return pd.DataFrame(rows)