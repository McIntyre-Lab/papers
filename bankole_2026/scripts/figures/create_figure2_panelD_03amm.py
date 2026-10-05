#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['savefig.dpi'] = 600

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY = f"{PROJ}/zenodo/summary_files"
ZENODO = f"{PROJ}/zenodo"
OUT = f"{PROJ}/Figures_amm"
os.makedirs(OUT, exist_ok=True)

species_list = ['dmel', 'dsim', 'dyak', 'dsan', 'dser']
species_to_genome = {'dmel': 'dmel6', 'dsim': 'dsim2', 'dyak': 'dyak2', 'dsan': 'dsan1', 'dser': 'dser1'}

species_colors_analyzable = {
    'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A', 'dyak': '#717273', 'dser': '#825CA6'
}
species_colors_unanalyzable = {
    'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7', 'dyak': '#B8B9BA', 'dser': '#C3A9D3'
}

transcript_flag_dict = {
    'monoExon': 'flag_transcriptClass_monoExon',
    'single_transcript': 'flag_transcriptClass_single',
    'multiple_transcripts': 'flag_transcriptClass_multiple'
}

comp_merged = pd.read_csv(f"{SUMMARY}/component_summary.csv", low_memory=False)
species_data = {}

for species in species_list:
    genome = species_to_genome[species]

    df_gene = pd.read_csv(
        f"{SUMMARY}/gene_summary_{species}.csv",
        usecols=['geneID', 'num_ujc', 'num_MonoExon', 'num_ujc_analyzable', 'num_ERPp_analyzable'], low_memory=False
    )

    df_gene['flag_analyzable'] = ((df_gene['num_ujc_analyzable'] > 0) | (df_gene['num_ERPp_analyzable'] > 0)).astype(int)
    if species in ['dyak', 'dsan']:
        df_gene['flag_analyzable'] = 0

    # monoExon: every annotated UJC has exactly 1 exon AND every observed UJC is mono-exon.
    anno = pd.read_csv(f"{ZENODO}/fiveSpecies_{genome}_full_annotation.csv", usecols=['geneID', 'numExon_ERP'], low_memory=False)
    anno['numExon_ERP'] = pd.to_numeric(anno['numExon_ERP'], errors='coerce')
    anno_all_mono = anno.groupby('geneID')['numExon_ERP'].apply(lambda x: x.notna().all() and (x == 1).all())
    df_gene['flag_all_anno_UJC_monoExon'] = df_gene['geneID'].map(anno_all_mono).fillna(False)
    
    #print((df_gene["num_ujc"] == 0).sum())
    
    df_gene['flag_all_data_UJC_monoExon'] = df_gene['num_ujc'] == df_gene['num_MonoExon']
    df_gene['flag_monoExon'] = df_gene['flag_all_anno_UJC_monoExon'] & df_gene['flag_all_data_UJC_monoExon']


    # single_transcript: if not monoExon, exactly one componentID and that component has flag_simple == 1.
    gene_col = f"geneID_concat_{species}"
    simple_col = f"flag_simple_{species}"
    comp_species = comp_merged[['componentID', gene_col, simple_col]].copy()
    comp_species['geneID'] = comp_species[gene_col].str.split('|')
    comp_expanded = comp_species.explode('geneID').dropna(subset=['geneID'])
    comp_expanded['geneID'] = comp_expanded['geneID'].str.strip()

    comp_gene = comp_expanded.groupby('geneID').agg(
        num_components=('componentID', 'nunique'),
        num_components_simple=(simple_col, 'sum')
    ).reset_index()
    comp_gene['flag_single_simple_component'] = (
        (comp_gene['num_components'] == 1) & (comp_gene['num_components_simple'] == 1)
    )

    single_simple_dict = comp_gene.set_index('geneID')['flag_single_simple_component'].to_dict()
    df_gene['flag_single_simple_component'] = df_gene['geneID'].map(single_simple_dict).fillna(False)

    df_gene['transcriptClass'] = np.select(
        [df_gene['flag_monoExon'], (~df_gene['flag_monoExon']) & df_gene['flag_single_simple_component']],
        ['monoExon', 'single_transcript'], default='multiple_transcripts'
    )

    df_gene['flag_transcriptClass_monoExon'] = (df_gene['transcriptClass'] == 'monoExon').astype(int)
    df_gene['flag_transcriptClass_single'] = (df_gene['transcriptClass'] == 'single_transcript').astype(int)
    df_gene['flag_transcriptClass_multiple'] = (df_gene['transcriptClass'] == 'multiple_transcripts').astype(int)

    species_data[species] = df_gene

category_order = list(transcript_flag_dict.keys())
count_data = {}

for species, df_gene in species_data.items():
    analyzable_conditions = {'analyzable': df_gene['flag_analyzable'] == 1, 'unanalyzable': df_gene['flag_analyzable'] == 0}
    transcript_conditions = {name: df_gene[flag] == 1 for name, flag in transcript_flag_dict.items()}
    count_data[species] = {}

    for analyzable_name, analyzable_condition in analyzable_conditions.items():
        count_data[species][analyzable_name] = {}
        for transcript_name, transcript_condition in transcript_conditions.items():
            count_data[species][analyzable_name][transcript_name] = int((analyzable_condition & transcript_condition).sum())


# Print counts as a species-by-count table
count_rows = {}
for transcript_name in category_order:
    count_rows[f"{transcript_name}_analyzable"] = [count_data[species]['analyzable'][transcript_name] for species in species_list]
    count_rows[f"{transcript_name}_unanalyzable"] = [count_data[species]['unanalyzable'][transcript_name] for species in species_list]
    count_rows[f"{transcript_name}_total"] = [count_data[species]['analyzable'][transcript_name] + count_data[species]['unanalyzable'][transcript_name] for species in species_list]

count_table = pd.DataFrame.from_dict(count_rows, orient='index', columns=['dmel', 'dsim', 'dyak', 'dsan', 'dser'])
count_table.index.name = 'counts'
print(count_table.to_csv(sep='\t'))

# count genes per species into plot
print("\nGENES PER SPECIES GOING INTO PLOT")
for species in species_list:
    total = sum(
        count_data[species]['analyzable'][cat] + count_data[species]['unanalyzable'][cat]
        for cat in category_order
    )
    print(f"  {species}: {total:,}")
    
    
plt.figure(figsize=(14, 7))
x = np.arange(len(category_order))
bar_width = 0.15
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for i, (species_name, counts) in enumerate(count_data.items()):
    x_pos = x + offset_positions[i] * bar_width
    analyzable_vals = [counts['analyzable'][cat] for cat in category_order]
    unanalyzable_vals = [counts['unanalyzable'][cat] for cat in category_order]
    plt.bar(x_pos, analyzable_vals, bar_width, color=species_colors_analyzable[species_name], edgecolor='black', linewidth=0.5)
    plt.bar(x_pos, unanalyzable_vals, bar_width, bottom=analyzable_vals,
            color=species_colors_unanalyzable[species_name], edgecolor='black', linewidth=0.5)

plt.xlabel('Transcript class', fontsize=12)
plt.ylabel('Number of Genes', fontsize=12)
plt.title('Genes by transcript class: analyzable vs unanalyzable', fontsize=14)
plt.xticks(x, category_order)
plt.ylim(0, 12000)  ## amm added y limit
plt.grid(axis='y', alpha=0.3, linestyle='--')
plt.tight_layout()

outfile = f"{OUT}/figure2_panelD_02mdg.pdf"
plt.savefig(outfile, format='pdf', bbox_inches='tight', dpi=600)
plt.close()
print(f"Saved: {outfile}")