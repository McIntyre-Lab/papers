#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 16:02:04 2025

@author: nkeil
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Function to read GTF file
def read_gtf(file_path, remove_decimals=False):
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    gtf_df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=gtf_columns)
    gtf_df['gene_id'] = gtf_df['attribute'].str.extract(r'gene_id "([^"]+)"')
    gtf_df['transcript_id'] = gtf_df['attribute'].str.extract(r'transcript_id "([^"]+)"')
    
    if remove_decimals:
        gtf_df['gene_id'] = gtf_df['gene_id'].apply(lambda x: str(x).split('.')[0] if pd.notnull(x) else x)
        gtf_df['transcript_id'] = gtf_df['transcript_id'].apply(lambda x: str(x).split('.')[0] if pd.notnull(x) else x)

    return gtf_df

# Function to read BED file
def read_bed(file_path):
    bed_columns = ['chrom', 'start', 'end', 'featureID']
    bed_df = pd.read_csv(file_path, sep='\t', header=None, names=bed_columns)
    bed_df['gene_id'] = bed_df['featureID'].apply(lambda x: x.split(':')[0])
    bed_df['transcript_id'] = 1  # Since only one transcript
    return bed_df

def plot_exons_by_transcript(df, gene_id, ax, color, exon_height=0.6, label=None, draw_connecting_lines=False, is_gtf=False, xlim=None):
    if is_gtf:
        gene_df = df[(df['gene_id'] == gene_id) & (df['feature'] == 'exon')]
    else:
        gene_df = df[df['gene_id'] == gene_id]
    
    # Drop duplicate exons
    gene_df = gene_df.drop_duplicates(subset=['start', 'end'])
    
    if gene_df.empty:
        return
    
    y_pos = 0
    prev_end = None
    prev_strand = None
    for _, row in gene_df.iterrows():
        ax.add_patch(Rectangle(
            (row['start'], y_pos - exon_height / 2), 
            row['end'] - row['start'], 
            exon_height, 
            facecolor=color, 
            edgecolor=color,
            alpha=1
        ))
        prev_end = row['end']
        prev_strand = row['strand'] if 'strand' in row else None
    
    ax.set_ylim(-exon_height / 2, exon_height / 2)
    if xlim:
        ax.set_xlim(xlim)
    ax.set_yticks([])

    if label:
        ax.annotate(label, xy=(0.5, 1), xytext=(0, 5), textcoords='offset points', ha='center', va='bottom')

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.xaxis.set_visible(False)

# File paths
mane_gtf_file_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/MANE.GRCh38.v1.2.ensembl_genomic.gtf'
bed_file_paths = [
    '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.bed',
    '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/hg38_union_event_analysis_ef.bed'
]
cd4_flag_file_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/flag_analyze_frag_DE_CD4.csv'
cd8_flag_file_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/flag_analyze_frag_DE_CD8.csv'

# Read the genes to be plotted
genes_file_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/t1d_gene_list_4_plots.csv'
genes_df = pd.read_csv(genes_file_path, header=None, names=['gene_name', 'gene_id'])

# Read BED files
all_frag_bed_dfs = [read_bed(path) for path in bed_file_paths]

# Read flag files
cd4_flag_df = pd.read_csv(cd4_flag_file_path)
cd8_flag_df = pd.read_csv(cd8_flag_file_path)

# Merge BED files with flag files
cd4_bed_df = pd.merge(all_frag_bed_dfs[0], cd4_flag_df, on='featureID')
cd8_bed_df = pd.merge(all_frag_bed_dfs[0], cd8_flag_df, on='featureID')

# Keep only rows with '1' in the last column
cd4_bed_df = cd4_bed_df[cd4_bed_df.iloc[:, -1] == 1]
cd8_bed_df = cd8_bed_df[cd8_bed_df.iloc[:, -1] == 1]

# Read GTF file
mane_gtf_df = read_gtf(mane_gtf_file_path, remove_decimals=True)

# Drop NaN values for transcript IDs in GTF file
mane_gtf_df = mane_gtf_df.dropna(subset=['transcript_id'])

# Process each gene
for index, row in genes_df.iterrows():
    gene_name = row['gene_name']
    gene_id = row['gene_id']
    
    # Filter all_frag BED file for exons longer than 10 nucleotides
    all_frag_bed_df = all_frag_bed_dfs[0][all_frag_bed_dfs[0]['gene_id'] == gene_id]
    all_frag_bed_df = all_frag_bed_df[all_frag_bed_df['end'] - all_frag_bed_dfs[0]['start'] > 10]  # <-- ensure >10nt filter

    # Filter merged BED files for the gene of interest
    cd4_bed_df_filtered = cd4_bed_df[cd4_bed_df['gene_id'] == gene_id]
    cd8_bed_df_filtered = cd8_bed_df[cd8_bed_df['gene_id'] == gene_id]

    # Filter the hg38 BED file for the gene of interest
    hg38_bed_df = all_frag_bed_dfs[1][all_frag_bed_dfs[1]['gene_id'] == gene_id]
    
    # Filter MANE GTF file for the gene of interest
    mane_gtf_df_filtered = mane_gtf_df[(mane_gtf_df['gene_id'] == gene_id) & (mane_gtf_df['feature'] == 'exon')]

    # Calculate the number of exons for each file
    num_exons_all_frag = len(all_frag_bed_df)
    num_exons_cd4 = len(cd4_bed_df_filtered)
    num_exons_cd8 = len(cd8_bed_df_filtered)
    num_exons_hg38 = len(hg38_bed_df)
    num_exons_mane = len(mane_gtf_df_filtered)

    # Debug: Print number of exons after filtering
    print(f'Number of exons in all_frag BED file for {gene_name} after filtering: {num_exons_all_frag}')
    print(f'Number of exons in CD4 BED file for {gene_name} after filtering: {num_exons_cd4}')
    print(f'Number of exons in CD8 BED file for {gene_name} after filtering: {num_exons_cd8}')
    print(f'Number of exons in hg38 BED file for {gene_name} after filtering: {num_exons_hg38}')
    print(f'Number of exons in MANE GTF file for {gene_name} after filtering: {num_exons_mane}')

    # Skip plotting if all DataFrames are empty
    if all(df.empty for df in [all_frag_bed_df, cd4_bed_df_filtered, cd8_bed_df_filtered, hg38_bed_df, mane_gtf_df_filtered]):
        print(f'No exons to plot for {gene_name}. Skipping...')
        continue

    # Determine x-axis limits based on all exons from all files
    starts = []
    ends = []
    for d in [all_frag_bed_df, cd4_bed_df_filtered, cd8_bed_df_filtered, hg38_bed_df, mane_gtf_df_filtered]:
        if not d.empty:
            starts.append(d['start'].min())
            ends.append(d['end'].max())
    if not starts or not ends:
        print(f'No valid exon bounds for {gene_name}. Skipping...')
        continue
    min_start = min(starts)
    max_end = max(ends)
    xlim = (min_start - 1000, max_end + 1000)

    # Plotting
    fig, axes = plt.subplots(5, 1, figsize=(20, 5), sharex=True)
    fig.suptitle(f'{gene_name} ({gene_id})')

    # Plot MANE GTF file
    plot_exons_by_transcript(mane_gtf_df_filtered, gene_id, axes[0], 'blue', label='MANE', is_gtf=True, xlim=xlim)

    # Plot hg38 BED file
    plot_exons_by_transcript(hg38_bed_df, gene_id, axes[1], 'orange', label='hg38_union', draw_connecting_lines=False, xlim=xlim)

    # Plot all_frag BED file
    plot_exons_by_transcript(all_frag_bed_df, gene_id, axes[2], 'red', label='All Frag', draw_connecting_lines=False, xlim=xlim)
    
    # Plot CD4 BED file
    plot_exons_by_transcript(cd4_bed_df_filtered, gene_id, axes[3], 'purple', label='CD4', draw_connecting_lines=False, xlim=xlim)

    # Plot CD8 BED file
    plot_exons_by_transcript(cd8_bed_df_filtered, gene_id, axes[4], 'green', label='CD8', draw_connecting_lines=False, xlim=xlim)

    # Make only the bottom plot show the x-axis
    axes[4].xaxis.set_visible(True)
    axes[4].spines['bottom'].set_visible(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for the title
    #plt.show()

    # Save the figure
    output_path = f'/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/transcript_plots_solid/{gene_name}.png'
    plt.savefig(output_path)

    plt.close(fig)
