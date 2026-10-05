#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['savefig.dpi'] = 600

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY = f"{PROJ}/zenodo/summary_files"
OUT = f"{PROJ}/Figures_amm"
os.makedirs(OUT, exist_ok=True)

species_list = ['dmel', 'dsim', 'dyak', 'dsan', 'dser']

species_colors_analyzable = {
    'dmel': '#966729', 'dsim': '#3F78C1', 'dsan': '#28827A',
    'dyak': '#717273', 'dser': '#825CA6'
}

species_colors_unanalyzable = {
    'dmel': '#D4B08C', 'dsim': '#9FC1E8', 'dsan': '#8FCCC7',
    'dyak': '#B8B9BA', 'dser': '#C3A9D3'
}

exon_flag_dict = {
    '1': 'flag_numExon_1',
    '2': 'flag_numExon_2',
    '3': 'flag_numExon_3',
    '4': 'flag_numExon_4',
    '5+': 'flag_numExon_5plus'
}

species_data = {}
for species in species_list:
    df_gene = pd.read_csv(
        f"{SUMMARY}/gene_summary_{species}.csv",
        usecols=['geneID', 'num_ujc_analyzable', 'num_ERPp_analyzable', 'numExon_GM'],
        low_memory=False
    )
    # flag analyzable if ujc OR ERPp is analyzed
    df_gene['flag_analyzable'] = (
        (df_gene['num_ujc_analyzable'] > 0) | (df_gene['num_ERPp_analyzable'] > 0)
    ).astype(int)

    # dyak and dsan were not evaluated for significant sex bias
    if species in ['dyak', 'dsan']:
        df_gene['flag_analyzable'] = 0
    
    # mutually exclusive categories based on numExon_GM 
    df_gene['flag_numExon_1'] = (df_gene['numExon_GM'] == 1).astype(int)
    df_gene['flag_numExon_2'] = (df_gene['numExon_GM'] == 2).astype(int)
    df_gene['flag_numExon_3'] = (df_gene['numExon_GM'] == 3).astype(int)
    df_gene['flag_numExon_4'] = (df_gene['numExon_GM'] == 4).astype(int)
    df_gene['flag_numExon_5plus'] = (df_gene['numExon_GM'] >= 5).astype(int)

    if df_gene["numExon_GM"].isna().any():
        print(f"{species}: has missing numExon_GM")
        ## check none missing

    ## check every gene falls in exactly one exon bin
    exon_flag_cols = list(exon_flag_dict.values())
    bin_sum = df_gene[exon_flag_cols].sum(axis=1)
    if not (bin_sum == 1).all():
        print(f"{species}: {(bin_sum == 0).sum()} genes in no bin, {(bin_sum > 1).sum()} genes in multiple bins")

    species_data[species] = df_gene

category_order = list(exon_flag_dict.keys())
count_data = {}

for species, df_gene in species_data.items():
    analyzable_conditions = {
        'analyzable': df_gene['flag_analyzable'] == 1,
        'unanalyzable': df_gene['flag_analyzable'] == 0
    }
    ## conditions for each 
    exon_conditions = {
        exon_name: df_gene[exon_flag] == 1
        for exon_name, exon_flag in exon_flag_dict.items()
    }

    count_data[species] = {}
    ## counts genes that are analyzable/unanalyzable for each exon category
    for analyzable_name, analyzable_condition in analyzable_conditions.items():
        count_data[species][analyzable_name] = {}
        for exon_name, exon_condition in exon_conditions.items():
            condition = analyzable_condition & exon_condition
            count_data[species][analyzable_name][exon_name] = int(condition.sum())

## print # genes per species to be plotted
for species in species_list:
    total = sum(
        count_data[species]['analyzable'][cat] + count_data[species]['unanalyzable'][cat]
        for cat in category_order
    )
    print(f"{species}: {total} genes plotted ")
        
    ## check plotted total matches original row count
    if total != len(species_data[species]):
        print(f"{species}: {len(species_data[species])} rows in df but {total} genes plotted")
#    print(f"{species}: {total} genes")


# Print counts as a species-by-count table
count_rows = {}
for exon_name in category_order:
    count_rows[f"{exon_name}_analyzable"] = [count_data[species]['analyzable'][exon_name] for species in species_list]
    count_rows[f"{exon_name}_unanalyzable"] = [count_data[species]['unanalyzable'][exon_name] for species in species_list]
    count_rows[f"{exon_name}_total"] = [count_data[species]['analyzable'][exon_name] + count_data[species]['unanalyzable'][exon_name] for species in species_list]

count_table = pd.DataFrame.from_dict(count_rows, orient='index', columns=['dmel', 'dsim', 'dyak', 'dsan', 'dser'])
count_table.index.name = 'counts'
print(count_table.to_csv(sep='\t'))

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

    
plt.xlabel('Number of Exons in gene model', fontsize=12)
plt.ylabel('Number of Genes', fontsize=12)
plt.title('Genes by exon number category: analyzable vs unanalyzable', fontsize=14)
plt.xticks(x, category_order)
plt.ylim(0, 12000)  ## amm changed to 12000 to match D
plt.grid(axis='y', alpha=0.3, linestyle='--')
plt.tight_layout()

outfile = f"{OUT}/figure2_panelC_02mdg.pdf"
plt.savefig(outfile, format='pdf', bbox_inches='tight', dpi=600)
plt.close()
print(f"Saved: {outfile}")