#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:41:39 2024

@author: nkeil
"""
import pandas as pd
import argparse
## Update options
def getOptions():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Make sqanti reads QC plots')
    
    parser.add_argument('-i','--in_file', type=str, dest="inFILE" ,required=True, help='Path to the T1D summary file')
    parser.add_argument('-o', '--out_path', type=str, dest="outPATH" ,required=True, help='Path to output')
    parser.add_argument('-p','--prefix', type=str, dest="outPREFIX" ,required=True, help='Output file prefix')
    
    
    args = parser.parse_args()
    return args


# Function to count and calculate proportions
def calc_counts_and_props(df, flag_sig_het_col, flag_sig_mod_col, immune_gene_col, t1d_gene_col):
    results = []

    categories = {
        'all_genes': df,
        'immune_genes': df[df[immune_gene_col] == 1],
        't1d_genes': df[df[t1d_gene_col] == 1]
    }

    for category_name, category_df in categories.items():
        total_genes = len(category_df)
        
        flag_sig_het_only = category_df[(category_df[flag_sig_het_col] == 1) & (category_df[flag_sig_mod_col] == 0)]
        flag_sig_mod_only = category_df[(category_df[flag_sig_het_col] == 0) & (category_df[flag_sig_mod_col] == 1)]
        flag_both = category_df[(category_df[flag_sig_het_col] == 1) & (category_df[flag_sig_mod_col] == 1)]
        
        counts = {
            'category': category_name,
            'analyzable_genes': total_genes,
            'splic_eff_only': len(flag_sig_het_only),
            'sex_eff_only': len(flag_sig_mod_only),
            'splic_and_sex_eff': len(flag_both)
        }

        props = {f'prop_{key}': value / total_genes if total_genes > 0 else 0 for key, value in counts.items() if key not in ['category', 'analyzable_genes']}
        
        results.append({**counts, **props})
    
    return pd.DataFrame(results)

def main():
    # Load the data
    file_path = args.inFILE  # Replace with your actual file path
    df = pd.read_csv(file_path)
    
    # Specify columns for each model
    models = {
    'overall': ('flag_sig_het_overall', 'flag_sig_mod_overall', 'flag_immune_gene', 'flag_t1d_gene'),
    'exons_only': ('flag_sig_het_exon', 'flag_sig_mod_exon', 'flag_immune_gene', 'flag_t1d_gene'),
    'introns_only': ('flag_sig_het_intron', 'flag_sig_mod_intron', 'flag_immune_gene', 'flag_t1d_gene')
}

    # Calculate counts and proportions for each model and store as DataFrames
    results = {}
    for model_name, columns in models.items():
        results[model_name] = calc_counts_and_props(df, *columns)
        
    out_path = args.outPATH
    prefix = args.outPREFIX

    for model_name, result in results.items():
        result.to_csv(f'{out_path}/{prefix}_{model_name}_model_counts.csv', index=False)

if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()