#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:04:44 2024

@author: nkeil
"""

import os
import re
import pandas as pd
import argparse

def getOptions():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process meta-analysis summary files and extract p-values")
    parser.add_argument("-d", "--directory", type=str, help="Path to the directory containing the files")
    parser.add_argument("-s", "--suffix", type=str, help="File suffix to filter files in the directory")
    parser.add_argument("-x", "--sig", type=float, help="P-value cutoff for significance")
    parser.add_argument("-o", "--out", type=str, help="Output CSV file name")
    
    args = parser.parse_args()
    return args

def extract_p_values(summary):
    # Flexible regex for heterogeneity and moderator p-values
    heterogeneity_pattern = r"Test for Residual Heterogeneity:.*?p-val\s*(<|=)\s*([.\d]+)"
    moderator_pattern = r"Test of Moderators.*?p-val\s*(<|=)\s*([.\d]+)"

    heterogeneity_match = re.search(heterogeneity_pattern, summary, re.DOTALL)
    moderator_match = re.search(moderator_pattern, summary, re.DOTALL)

    def format_p_value(match):
        if match:
            operator, value = match.groups()
            return 0.0001 if operator == "<" else float(value)
        else:
            return None

    p_val_heterogeneity = format_p_value(heterogeneity_match)
    p_val_moderator = format_p_value(moderator_match)

    return p_val_heterogeneity, p_val_moderator

def process_files(directory, suffix):
    results = []

    for filename in os.listdir(directory):
        if filename.endswith(suffix):
            print(f"Processing file: {filename}")
            with open(os.path.join(directory, filename), 'r') as file:
                summary = file.read()

            p_heterogeneity, p_moderator = extract_p_values(summary)
            if (p_heterogeneity is not None) and (p_moderator is not None):
                filename_base = os.path.basename(filename)
                parts = filename_base.split("_")
                
                # Search for CD4 or CD8 dynamically
                if "CD4" in parts or "CD8" in parts:
                    try:
                        idx = parts.index("CD4") if "CD4" in parts else parts.index("CD8")
                        gene_id = parts[idx-1]  # Gene ID is the part just before CD4/CD8
                        cell_type = parts[idx]

                        results.append([gene_id, cell_type, p_heterogeneity, p_moderator])
                    except IndexError:
                        print(f"⚠️ Warning: Unexpected filename structure: {filename}")
                else:
                    print(f"⚠️ Warning: No cell type CD4/CD8 found in {filename}")
            else:
                print(f"⚠️ No p-values found in {filename}")

    return pd.DataFrame(results, columns=['geneID', 'cellType', 'pval_heterogeneity', 'pval_moderator'])

def main():
    pvalDf = process_files(args.directory, args.suffix)

    pvalDf['flag_sig_heterogeneity'] = (pvalDf['pval_heterogeneity'] < args.sig).astype(int)
    pvalDf['flag_sig_moderator'] = (pvalDf['pval_moderator'] < args.sig).astype(int)
    
    pvalDf.to_csv(args.out, index=False)
    print(f"✅ Output saved to {args.out}")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()