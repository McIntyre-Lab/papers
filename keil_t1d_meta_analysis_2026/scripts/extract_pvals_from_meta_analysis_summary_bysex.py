#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 15:11:42 2025

@author: nkeil
"""
import os
import re
import pandas as pd
import argparse

def getOptions():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Process metanalysis summary files and extract p-values")
    parser.add_argument("-d", "--directory", type=str, help="Path to the directory containing the files")
    parser.add_argument("-s", "--suffix", type=str, help="File suffix to filter files in the directory")
    parser.add_argument("-x", "--sig", type=float, help="Pvalue cutoff for significance")
    parser.add_argument("-o", "--out", type=str, help="Output CSV file name")
    
    args = parser.parse_args()
    return args

def extract_p_values(summary):
    # Updated regex: captures both "<" or "=" and the number (even .0001)
    heterogeneity_pattern = r"Test for Heterogeneity:.*?p-val\s*(<|=)\s*([.\d]+)"
    
    heterogeneity_match = re.search(heterogeneity_pattern, summary, re.DOTALL)

    def format_p_value(match):
        if match:
            operator, value = match.groups()
            return 0.0001 if operator == "<" else float(value)
        else:
            return None

    return format_p_value(heterogeneity_match)

def process_files(directory, suffix):
    results = []

    for filename in os.listdir(directory):
        if filename.endswith(suffix):
            print(f"Processing file: {filename}")
            with open(os.path.join(directory, filename), 'r') as file:
                summary = file.read()

            p_heterogeneity = extract_p_values(summary)
            if p_heterogeneity is not None:
                filename_base = os.path.basename(filename)
                parts = filename_base.split("_")
                
                if len(parts) >= 4:
                    # Assume standard structure: *_GENEID_CELLTYPE_SEX_*
                    gene_id = parts[1]  # GENEID (e.g., ENSG00000001630)
                    cell_type = parts[2]  # CELLTYPE (e.g., CD4)
                    sex = parts[3]  # SEX (e.g., F, M)
                    
                    # Safety check: is cell_type actually CD4 or CD8?
                    if cell_type not in ["CD4", "CD8"]:
                        print(f"⚠️ Warning: Unexpected cell type in {filename}: got {cell_type}")
                        cell_type = "Unknown"

                    # Safety check: is sex either F or M?
                    if sex not in ["F", "M"]:
                        print(f"⚠️ Warning: Unexpected sex in {filename}: got {sex}")
                        sex = "Unknown"

                    results.append([gene_id, cell_type, sex, p_heterogeneity])
                else:
                    print(f"⚠️ Warning: Unexpected filename format: {filename}")
            else:
                print(f"⚠️ No p-values found in {filename}")

    return pd.DataFrame(results, columns=['geneID', 'cellType', 'sex', 'pval_heterogeneity'])

def main():
    # Get options from command line args
    pvalDf = process_files(args.directory, args.suffix)
    
    pvalDf['flag_sig_heterogeneity'] = (pvalDf['pval_heterogeneity'] < args.sig).astype(int)
    
    # Save the results to a CSV file
    pvalDf.to_csv(args.out, index=False)
    print(f"Output saved to {args.out}")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
