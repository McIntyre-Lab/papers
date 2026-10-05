#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_component_network_graph_02mdg.py

Rendering module for component network graphs.

The core function plot_component_network() accepts pre-filtered DataFrames
and handles only graph construction and drawing — no file I/O.

The CLI accepts the full node map and edges files and filters to the
requested component_id internally, so the bash caller does not need to
pre-filter anything.

Node coloring
    effect_sizes provided : blue-purple-red gradient, species label on node
    effect_sizes absent   : SPECIES_COLORS by genome string, species legend

Edge coloring
    Single uniform color supplied by the caller (e.g. red / blue / gray).

Layout
    Identical to plot_component_network_graph_01mdg.py:
    core nodes (degree >= 2) in a circle,
    leaf nodes (degree == 1) radially outward from their parent.
    Graph type: nx.MultiDiGraph.
"""

import argparse
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SPECIES_COLORS = {
    'dmel6': '#966729',
    'dsim2': '#3F78C1',
    'dsan1': '#28827A',
    'dyak2': '#717273',
    'dser1': '#825CA6',
}

SPECIES_SHORT = {
    'dmel6': 'mel',
    'dsim2': 'sim',
    'dsan1': 'san',
    'dyak2': 'yak',
    'dser1': 'ser',
}

EFFECT_CMAP = LinearSegmentedColormap.from_list(
    'effect_bpr',
    [(0.0, (0.0, 0.0, 1.0)),
     (0.5, (0.5, 0.0, 0.5)),
     (1.0, (1.0, 0.0, 0.0))],
)
EFFECT_VMIN   = -8
EFFECT_VMAX   =  8
NO_DATA_COLOR = (0.75, 0.75, 0.75)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def effect_to_color(effect_size):
    if pd.isna(effect_size):
        return NO_DATA_COLOR
    norm_val = np.clip(
        (float(effect_size) - EFFECT_VMIN) / (EFFECT_VMAX - EFFECT_VMIN),
        0.0, 1.0,
    )
    return EFFECT_CMAP(norm_val)


def build_layout(G):
    """Geometric layout: core nodes in a circle, leaves radially outward."""
    core_nodes = [n for n in G.nodes() if G.degree(n) >= 2]
    leaf_nodes = [n for n in G.nodes() if G.degree(n) == 1]

    pos = {}
    n_core = len(core_nodes)
    for i, node in enumerate(core_nodes):
        angle = 2 * np.pi * i / n_core if n_core > 1 else 0
        pos[node] = (np.cos(angle), np.sin(angle))

    for leaf in leaf_nodes:
        neighbors = list(G.neighbors(leaf)) + list(G.predecessors(leaf))
        if not neighbors:
            pos[leaf] = (0, 0)
            continue
        parent = neighbors[0]
        if parent in pos:
            parent_pos = np.array(pos[parent])
            norm      = np.linalg.norm(parent_pos)
            direction = parent_pos / norm if norm > 0 else np.array([1, 0])
            pos[leaf] = tuple(parent_pos + direction * 0.5)
        else:
            pos[leaf] = (0, 0)

    return pos


# ---------------------------------------------------------------------------
# Core plotting function
# ---------------------------------------------------------------------------

def plot_component_network(
    component_id,
    nodes,
    edges,
    output_dir,
    effect_sizes=None,
    edge_color='#333333',
    title=None,
    output_format='png',
):
    """
    Render and save one network graph for a single component.

    Parameters
    ----------
    component_id  : int or str
    nodes         : DataFrame with columns [jxnhash, source]
                    Already filtered to this component.
    edges         : DataFrame with columns [source_jxnHash, target_jxnHash]
                    Already filtered to this component.
    output_dir    : str
    effect_sizes  : dict  jxnhash -> float, or None
    edge_color    : str, any matplotlib color
    title         : str or None
    output_format : 'png' | 'svg' | 'pdf'
    """
    use_effect = effect_sizes is not None

    G = nx.MultiDiGraph()
    for _, row in nodes.iterrows():
        es = effect_sizes.get(row['jxnhash'], np.nan) if use_effect else None
        G.add_node(
            row['jxnhash'],
            species       = row['source'],
            effect_size   = es,
            species_label = SPECIES_SHORT.get(row['source'], row['source'][:3]),
        )

    node_set = set(G.nodes())
    for _, row in edges.iterrows():
        src, tgt = row['source_jxnHash'], row['target_jxnHash']
        if src in node_set and tgt in node_set:
            G.add_edge(src, tgt)

    pos       = build_layout(G)
    node_list = list(G.nodes())

    if use_effect:
        node_colors = [effect_to_color(G.nodes[n]['effect_size']) for n in node_list]
        node_labels = {n: G.nodes[n]['species_label'] for n in node_list}
    else:
        node_colors = [SPECIES_COLORS.get(G.nodes[n]['species'], '#CCCCCC')
                       for n in node_list]
        node_labels = None

    fig, ax = plt.subplots(figsize=(6, 6))

    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color=edge_color,
                           arrows=True, arrowsize=15)
    nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=node_list,
                           node_color=node_colors, node_size=300, alpha=0.9)

    if use_effect:
        nx.draw_networkx_labels(G, pos, ax=ax, labels=node_labels,
                                font_size=6, font_color='white', font_weight='bold')
        sm = ScalarMappable(cmap=EFFECT_CMAP,
                            norm=Normalize(vmin=EFFECT_VMIN, vmax=EFFECT_VMAX))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.03, aspect=20)
        cbar.set_label('effect_size_Equal', fontsize=9)
        cbar.set_ticks([EFFECT_VMIN, 0, EFFECT_VMAX])
    else:
        species_counts  = nodes['source'].value_counts()
        legend_elements = [
            Patch(facecolor=SPECIES_COLORS[sp],
                  label=f'{sp} (n={species_counts[sp]})')
            for sp in SPECIES_COLORS
            if sp in species_counts.index
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=7)

    ax.set_title(title or f'Component {component_id}', fontsize=11, fontweight='bold')
    ax.axis('off')
    ax.set_aspect('equal')
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir,
                            f'component_{component_id}_network.{output_format}')
    save_kw = {} if output_format == 'svg' else {'dpi': 200}
    fig.savefig(out_path, bbox_inches='tight', **save_kw)
    plt.close(fig)
    print(f'saved: {out_path}')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Render one component network graph.'
    )
    parser.add_argument('--component_id',      required=True)
    parser.add_argument('--nodes_file',        required=True,
                        help='Full component_map_by_node.csv '
                             '(filtered to component_id internally)')
    parser.add_argument('--edges_file',        required=True,
                        help='Full edges.csv '
                             '(filtered to component_id internally)')
    parser.add_argument('--output_dir',        required=True)
    parser.add_argument('--effect_sizes_file', default=None,
                        help='CSV with columns: jxnhash, effect_size_Equal')
    parser.add_argument('--edge_color',        default='#333333')
    parser.add_argument('--title',             default=None)
    parser.add_argument('--output_format',     choices=['png', 'svg', 'pdf'],
                        default='png')
    args = parser.parse_args()

    nodes_df = pd.read_csv(args.nodes_file, dtype=str)
    edges_df = pd.read_csv(args.edges_file, dtype=str)

    # Filter to the requested component
    nodes_df = nodes_df[nodes_df['component_id'].astype(str) == str(args.component_id)]
    node_set = set(nodes_df['jxnhash'])
    edges_df = edges_df[
        edges_df['source_jxnHash'].isin(node_set) &
        edges_df['target_jxnHash'].isin(node_set)
    ]

    es = None
    if args.effect_sizes_file:
        es_df = pd.read_csv(args.effect_sizes_file)
        es    = dict(zip(es_df['jxnhash'], es_df['effect_size_Equal']))

    plot_component_network(
        component_id  = args.component_id,
        nodes         = nodes_df,
        edges         = edges_df,
        output_dir    = args.output_dir,
        effect_sizes  = es,
        edge_color    = args.edge_color,
        title         = args.title,
        output_format = args.output_format,
    )