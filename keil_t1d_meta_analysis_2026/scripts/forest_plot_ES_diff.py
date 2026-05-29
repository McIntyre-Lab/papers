#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 21 13:32:28 2025

@author: nkeil
"""
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# Path to the directory containing the CSV files
dir_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/data_4_forest_plots_ES_diff/'
graphs_path = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/test_files_for_anna/graphs/forest_plots_ES_diff'

# Create directory for graphs if it doesn't exist
os.makedirs(graphs_path, exist_ok=True)

# Plotting function
def plot_data(df, gene_name, cell_type, title, ylim, filename):
    num_fragments = len(df)
    plt.figure(figsize=(24, 8))  # Adjust width based on number of fragments
    x_positions = range(num_fragments)

    # Add sample plot elements for the legend
    plt.plot([], [], 's', color='black', alpha=1.0, label='Exon')
    plt.plot([], [], 's', color='black', alpha=0.5, linestyle='dashed', label='Intron')

    # Plot for all data
    for x, row in zip(x_positions, df.itertuples()):
        alpha = 0.5 if row.ef_ir_flag == 1 else 1.0
        linestyle = 'dashed' if row.ef_ir_flag == 1 else 'solid'
        # Error bar plot
        (line, caps, bars) = plt.errorbar(x, row.effect_size, 
                                          yerr=[[row.effect_size - row.ci_lower], [row.ci_upper - row.effect_size]],
                                          fmt='s', capsize=5, ecolor='black', color='black', markeredgecolor='black', alpha=alpha)
        for bar in bars:
            bar.set_linestyle(linestyle)  # Apply dashed linestyle if ef_ir_flag == 1
            bar.set_alpha(alpha)

    plt.axhline(0, color='grey', linestyle='--', alpha=0.5)  # Add grey transparent dotted line at y=0
    plt.xlabel('Exon Fragments')
    plt.ylabel('Effect Size Difference (Male - Female )', labelpad=20)  # Ensure y-axis title is fully displayed
    plt.title(f' Effect Size Difference (Male - Female) for {gene_name} ({cell_type}) by Fragment')
    plt.grid(True, alpha=0.5)
    plt.xlim(-0.5, num_fragments - 0.5)  # Ensure even spacing
    plt.ylim(ylim)
    plt.xticks(ticks=x_positions, labels=df['featureID'], rotation=90)  # Label x-axis with featureID
    plt.legend()
    plt.tight_layout(pad=3.0, w_pad=3.0, h_pad=3.0)  # Adjust layout to ensure the y-axis title is fully displayed
    plt.savefig(os.path.join(graphs_path, filename))  # Save the plot
    plt.show()
    plt.close()

# Loop through each file in the directory
for file in os.listdir(dir_path):
    if '_ES_diff_' and 'features' in file and file.endswith('.csv'):
        file_path = os.path.join(dir_path, file)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Extract the gene name and cell type from the filename
        gene_name = re.search(r'ENSG\d+', file).group(0)
        cell_type_match = re.search(r'CD[48]', file)
        cell_type = cell_type_match.group(0) if cell_type_match else 'Unknown'
        
        # Sort the DataFrame based on the start position
        df = df.sort_values(by='ef_start').reset_index(drop=True)
        
        # Create a new column 'featureID' excluding the gene name for the x-axis
        df['featureID'] = df['featureID'].apply(lambda x: ':'.join(x.split(':')[1:]))
        
        # Determine y-axis range based on the actual bounds
        min_y = df['ci_lower'].min()
        max_y = df['ci_upper'].max()
        y_range = (min_y - 0.25, max_y + 0.25)
        
        # Plot for all data
        plot_data(df, gene_name, cell_type, 'Effect Size by Fragment', ylim=y_range, filename=f'ef_sz_diff_{gene_name}_{cell_type}.png')
