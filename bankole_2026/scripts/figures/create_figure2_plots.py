#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import defaultdict

# Set global font to Arial for publication consistency
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'

BASE_DIR    = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures"
SUMMARY     = f"{BASE_DIR}/zenodo/summary_files"
NETWORK     = f"{BASE_DIR}/zenodo/fiveSpecies_supporting_files/network_files"
TABLES      = f"{BASE_DIR}/Tables"

###############################################################################
# PANEL A (Left): comp_anno_in value counts for expressed components

comp_df = pd.read_csv(f"{NETWORK}/component_info.csv")

category_map = {
    '4_yakonly':       'one_species_yak/mel/san/sim',
    '4_melonly':       'one_species_yak/mel/san/sim',
    '4_sanonly':       'one_species_yak/mel/san/sim',
    '4_simonly':       'one_species_yak/mel/san/sim',
    '4_zseronly':      'one_species_ser',
    '3_mel_sim':       'other',
    '3_yak_san':       'other',
    '6_other':         'other',
    '2_the4':          'mel,sim,san,yak',
    '1_all_5_species': 'mel,sim,san,yak,ser'
}

comp_df['flag_panelA_left'] = (comp_df['flag_expressed'] == 1).astype(int)
raw_cnts = comp_df.groupby('comp_anno_in')['flag_panelA_left'].sum()
mapped_index = pd.Index([category_map.get(x, x) for x in raw_cnts.index])
combined_cnts = pd.Series(index=mapped_index, data=raw_cnts.values).groupby(level=0).sum()

order = ['one_species_yak/mel/san/sim', 'one_species_ser', 'other',
         'mel,sim,san,yak', 'mel,sim,san,yak,ser']
combined_cnts = combined_cnts.reindex(order, fill_value=0)

plt.figure(figsize=(10, 8))
plt.pie(combined_cnts.values, labels=combined_cnts.index, autopct='%1.1f%%',
        startangle=210, counterclock=False,
        colors=['#E3F2FD', '#90CAF9', '#42A5F5', '#1E88E5', '#0D47A1'])
plt.title('Expressed components annotated in (comp_anno_in)')
plt.axis('equal')
plt.savefig(f"{FIGURES_DIR}/figure2_panelA_left.svg", format='svg', bbox_inches='tight')
plt.close()
print(f"Saved: {FIGURES_DIR}/figure2_panelA_left.svg")

###############################################################################
# PANEL A (Center): Sankey diagram (Comp Anno -> Expressed In)

sankey_df = pd.read_csv(f"{TABLES}/cnt_origin_anno_xpress_byspecies.csv")
sankey_df['flag_panelA_sankey'] = (
    sankey_df['expressed_in'] != '7_not_expressed'
).astype(int)

keep = ['1_all_5_species', '2_the4', '4_zseronly']
sankey_df['comp_anno_in'] = sankey_df['comp_anno_in'].where(
    sankey_df['comp_anno_in'].isin(keep), 'other'
)

comp_anno_flags = sankey_df.groupby('comp_anno_in')['flag_panelA_sankey'].max().to_dict()
expressed_flags = sankey_df.groupby('expressed_in')['flag_panelA_sankey'].max().to_dict()
comp_annos = [cat for cat, flag in comp_anno_flags.items() if flag == 1]
expressed_vals = [cat for cat, flag in expressed_flags.items() if flag == 1]

comp_anno_labels = [f"Comp Anno: {x}" for x in comp_annos]
expressed_labels = [f"Expressed: {x}" for x in expressed_vals]
all_labels = comp_anno_labels + expressed_labels
label_to_idx = {label: idx for idx, label in enumerate(all_labels)}

link_counts = defaultdict(int)
for _, row in sankey_df.iterrows():
    if row['flag_panelA_sankey'] == 1:
        source_label = f"Comp Anno: {row['comp_anno_in']}"
        target_label = f"Expressed: {row['expressed_in']}"
        source_idx = label_to_idx[source_label]
        target_idx = label_to_idx[target_label]
        link_counts[(source_idx, target_idx)] += row['COUNT']

sources = [k[0] for k in link_counts]
targets = [k[1] for k in link_counts]
values = list(link_counts.values())

fig_sankey = go.Figure(data=[go.Sankey(
    node=dict(pad=15, thickness=20, line=dict(color='black', width=0.5), label=all_labels, color='blue'),
    link=dict(source=sources, target=targets, value=values)
)])
fig_sankey.update_layout(title="Sankey Plot: Comp Anno → Expressed", font=dict(size=10), height=800, width=1200)
sankey_output_path = f"{FIGURES_DIR}/figure2_panelA_sankey_anno_2_expr.svg"
fig_sankey.write_image(sankey_output_path, format='svg')
print(f"Saved: {sankey_output_path}")

###############################################################################
# PANEL A (Right): expressed_in value counts for components annotated in all 5 species

species_colors = {
    'mel': '#966729', 'sim': '#3F78C1', 'san': '#28827A',
    'yak': '#717273', 'ser': '#825CA6'
}

comp_df['flag_panelA_right'] = (
    (comp_df['flag_expressed'] == 1)
    & (comp_df['comp_anno_in'] == '1_all_5_species')
).astype(int)
expressed_in_count_dict = comp_df.groupby('expressed_in')['flag_panelA_right'].sum().to_dict()

order_right = ['4_mel_only', '4_sim_only', '4_san_only', '4_yak_only', '4_zser_only',
               '3_mel_sim', '3_yak_san', '6_other', '2_the4', '1_all_5']
expressed_in_cnts = pd.Series({
    cat: expressed_in_count_dict.get(cat, 0)
    for cat in order_right
    if expressed_in_count_dict.get(cat, 0) > 0
})

style_map = {
    '4_mel_only':  (species_colors['mel'], 'white', ''),
    '4_sim_only':  (species_colors['sim'], 'white', ''),
    '4_san_only':  (species_colors['san'], 'white', ''),
    '4_yak_only':  (species_colors['yak'], 'white', ''),
    '4_zser_only': (species_colors['ser'], 'white', ''),
    '3_mel_sim':   (species_colors['mel'], species_colors['sim'], '//////'),
    '3_yak_san':   (species_colors['san'], species_colors['yak'], '//////'),
    '2_the4':      ('#1E88E5', 'white', ''),
    '6_other':     ('#42A5F5', 'white', ''),
    '1_all_5':     ('#0D47A1', 'white', ''),
}
colors     = [style_map.get(cat, ('#BBBBBB', 'white', ''))[0] for cat in expressed_in_cnts.index]
edgecolors = [style_map.get(cat, ('#BBBBBB', 'white', ''))[1] for cat in expressed_in_cnts.index]
hatches    = [style_map.get(cat, ('#BBBBBB', 'white', ''))[2] for cat in expressed_in_cnts.index]
labels     = [f"{label}\n({val / expressed_in_cnts.sum() * 100:.1f}%)"
              for label, val in zip(expressed_in_cnts.index, expressed_in_cnts.values)]

fig, ax = plt.subplots(figsize=(10, 8))
wedges, _ = ax.pie(expressed_in_cnts.values, labels=labels,
                    startangle=230, counterclock=False, colors=colors)
for wedge, hatch, ec in zip(wedges, hatches, edgecolors):
    wedge.set_hatch(hatch)
    wedge.set_edgecolor(ec)
    wedge.set_linewidth(0.0)
plt.title('Expressed components annotated in all 5 species by expression profile')
plt.axis('equal')
plt.savefig(f"{FIGURES_DIR}/figure2_panelA_right.svg", format='svg', bbox_inches='tight')
plt.close()
print(f"Saved: {FIGURES_DIR}/figure2_panelA_right.svg")

###############################################################################
# PANELS B AND C DATA PREP

species_list = ['dmel', 'dsim', 'dyak', 'dsan', 'dser']

species_colors_analyzable = {
    'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A',
    'dyak': '#717273', 'dser': '#825CA6'
}
species_colors_unanalyzable = {
    'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7',
    'dyak': '#B8B9BA', 'dser': '#C3A9D3'
}

comp_merged = pd.read_csv(f"{SUMMARY}/component_summary.csv", low_memory=False)

exon_flag_dict = {
    '1': 'flag_numExon_1',
    '2': 'flag_numExon_2',
    '3': 'flag_numExon_3',
    '4': 'flag_numExon_4',
    '5+': 'flag_numExon_5plus'
}

transcript_flag_dict = {
    'monoExon': 'flag_transcriptClass_monoExon',
    'single_transcript': 'flag_transcriptClass_single',
    'multiple_transcripts': 'flag_transcriptClass_multiple'
}

species_data = {}
for species in species_list:

    # Read the gene summary and flag genes that are analyzable
    df_gene = pd.read_csv(
        f"{SUMMARY}/gene_summary_{species}.csv",
        usecols=['geneID', 'num_ujc_analyzable', 'num_ERPp_analyzable', 'numExon_GM'],
        low_memory=False
    )

    df_gene['flag_analyzable'] = (
        (df_gene['num_ujc_analyzable'] > 0)
        | (df_gene['num_ERPp_analyzable'] > 0)
    ).astype(int)

    # Flag exon-number categories using numExon_GM from the gene summary
    df_gene['flag_numExon_1'] = (df_gene['numExon_GM'] == 1).astype(int)
    df_gene['flag_numExon_2'] = (df_gene['numExon_GM'] == 2).astype(int)
    df_gene['flag_numExon_3'] = (df_gene['numExon_GM'] == 3).astype(int)
    df_gene['flag_numExon_4'] = (df_gene['numExon_GM'] == 4).astype(int)
    df_gene['flag_numExon_5plus'] = (df_gene['numExon_GM'] >= 5).astype(int)

    # Expand species-specific concatenated gene IDs in the component summary
    gene_col = f"geneID_concat_{species}"
    simple_col = f"flag_simple_{species}"
    monoExon_col = f"flag_monoExon_{species}"

    comp_species = comp_merged[
        ['componentID', gene_col, simple_col, monoExon_col]
    ].copy()

    comp_species['geneID'] = comp_species[gene_col].str.split('|')
    comp_expanded = comp_species.explode('geneID')
    comp_expanded['geneID'] = comp_expanded['geneID'].str.strip()

    # Count total, simple, and monoExon components for each gene
    comp_gene = comp_expanded.groupby('geneID').agg(
        num_components=('componentID', 'nunique'),
        num_components_simple=(simple_col, 'sum'),
        num_components_monoExon=(monoExon_col, 'sum')
    ).reset_index()

    # Classify genes in priority order:
    # 1) all components monoExon
    # 2) otherwise exactly one simple component
    # 3) otherwise multiple transcripts
    flag_all_monoExon = (
        comp_gene['num_components_monoExon'] == comp_gene['num_components']
    )
    flag_one_simple = comp_gene['num_components_simple'] == 1

    comp_gene['transcriptClass'] = np.select(
        [flag_all_monoExon, flag_one_simple],
        ['monoExon', 'single_transcript'],
        default='multiple_transcripts'
    )

    # Add transcriptClass to the gene summary without merging or concatenating
    transcriptClass_dict = comp_gene.set_index('geneID')['transcriptClass'].to_dict()
    df_gene['transcriptClass'] = df_gene['geneID'].map(transcriptClass_dict)

    df_gene['flag_transcriptClass_monoExon'] = (
        df_gene['transcriptClass'] == 'monoExon'
    ).astype(int)
    df_gene['flag_transcriptClass_single'] = (
        df_gene['transcriptClass'] == 'single_transcript'
    ).astype(int)
    df_gene['flag_transcriptClass_multiple'] = (
        df_gene['transcriptClass'] == 'multiple_transcripts'
    ).astype(int)

    print(
        species,
        "genes:", len(df_gene),
        "analyzable:", int(df_gene['flag_analyzable'].sum()),
        "missing transcriptClass:", int(df_gene['transcriptClass'].isna().sum())
    )

    species_data[species] = df_gene

###############################################################################
# PANEL B: genes by number of exons in the gene model

category_order_b = list(exon_flag_dict.keys())
count_data_b = {}

for species, df_gene in species_data.items():

    analyzable_conditions = {
        'analyzable': df_gene['flag_analyzable'] == 1,
        'unanalyzable': df_gene['flag_analyzable'] == 0
    }
    exon_conditions = {
        exon_name: df_gene[exon_flag] == 1
        for exon_name, exon_flag in exon_flag_dict.items()
    }

    count_data_b[species] = {}
    for analyzable_name, analyzable_condition in analyzable_conditions.items():
        count_data_b[species][analyzable_name] = {}
        for exon_name, exon_condition in exon_conditions.items():
            condition = analyzable_condition & exon_condition
            count_data_b[species][analyzable_name][exon_name] = int(condition.sum())

plt.figure(figsize=(14, 7))
x = np.arange(len(category_order_b))
bar_width = 0.15
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for i, (species_name, counts) in enumerate(count_data_b.items()):
    x_pos = x + offset_positions[i] * bar_width
    analyzable_vals = [counts['analyzable'][cat] for cat in category_order_b]
    unanalyzable_vals = [counts['unanalyzable'][cat] for cat in category_order_b]

    plt.bar(
        x_pos, analyzable_vals, bar_width,
        color=species_colors_analyzable[species_name],
        edgecolor='black', linewidth=0.5
    )
    plt.bar(
        x_pos, unanalyzable_vals, bar_width,
        bottom=analyzable_vals,
        color=species_colors_unanalyzable[species_name],
        edgecolor='black', linewidth=0.5
    )

plt.xlabel('Number of Exons in gene model', fontsize=12)
plt.ylabel('Number of Genes', fontsize=12)
plt.title('Genes by exon number category: analyzable vs unanalyzable', fontsize=14)
plt.xticks(x, category_order_b)
plt.ylim(0, 10000)
plt.grid(axis='y', alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure2_panelB.svg", format='svg', bbox_inches='tight')
plt.close()
print(f"Saved: {FIGURES_DIR}/figure2_panelB.svg")

###############################################################################
# PANEL C: genes by transcriptClass

category_order_c = list(transcript_flag_dict.keys())
count_data_c = {}

for species, df_gene in species_data.items():

    analyzable_conditions = {
        'analyzable': df_gene['flag_analyzable'] == 1,
        'unanalyzable': df_gene['flag_analyzable'] == 0
    }
    transcript_conditions = {
        transcript_name: df_gene[transcript_flag] == 1
        for transcript_name, transcript_flag in transcript_flag_dict.items()
    }

    count_data_c[species] = {}
    for analyzable_name, analyzable_condition in analyzable_conditions.items():
        count_data_c[species][analyzable_name] = {}
        for transcript_name, transcript_condition in transcript_conditions.items():
            condition = analyzable_condition & transcript_condition
            count_data_c[species][analyzable_name][transcript_name] = int(condition.sum())

plt.figure(figsize=(14, 7))
x = np.arange(len(category_order_c))
bar_width = 0.15
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for i, (species_name, counts) in enumerate(count_data_c.items()):
    x_pos = x + offset_positions[i] * bar_width
    analyzable_vals = [counts['analyzable'][cat] for cat in category_order_c]
    unanalyzable_vals = [counts['unanalyzable'][cat] for cat in category_order_c]

    plt.bar(
        x_pos, analyzable_vals, bar_width,
        color=species_colors_analyzable[species_name],
        edgecolor='black', linewidth=0.5
    )
    plt.bar(
        x_pos, unanalyzable_vals, bar_width,
        bottom=analyzable_vals,
        color=species_colors_unanalyzable[species_name],
        edgecolor='black', linewidth=0.5
    )

plt.xlabel('Transcript class', fontsize=12)
plt.ylabel('Number of Genes', fontsize=12)
plt.title('Genes by transcript class: analyzable vs unanalyzable', fontsize=14)
plt.xticks(x, category_order_c)
plt.grid(axis='y', alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure2_panelC.svg", format='svg', bbox_inches='tight')
plt.close()
print(f"Saved: {FIGURES_DIR}/figure2_panelC.svg")