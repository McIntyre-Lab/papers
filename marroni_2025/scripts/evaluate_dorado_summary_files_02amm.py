#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 16:25:26 2024

@author: ammorse
"""
import os
import pandas as pd
import argparse

def process_summary_file(input_dir, output_dir, file_name, combined_match_counts):
    file_path = os.path.join(input_dir, file_name)
    df = pd.read_csv(file_path, sep="\t", header=0)

    df['input_pass_fail'] = df['filename'].apply(lambda x: x.split('_')[0])
    df['input_bc'] = df['filename'].apply(lambda x: x.split('_')[2])
    df['input_barcode'] = df['input_bc'].apply(lambda x: x.split('.')[0])

    df['output_barcode'] = df.apply(lambda row: row['barcode'].split('_')[1] if 
                                'SQK' in row['barcode'] else row['barcode'], 
                                axis=1)

    df.drop(columns=['input_bc'], inplace=True)

    ## check read_id uniq
    cnt_id_df = df['read_id'].value_counts().reset_index()
    cnt_id_df.columns = ['read_id', 'count']

    df['match'] = 'mismatch'  # Initialize the 'match' column with 'mismatch'

    # Update 'match' column based on conditions
    df.loc[df['input_barcode'] == df['output_barcode'], 'match'] = 'match'
    df.loc[df['input_barcode'] == 'unclassified', 'match'] = 'new_assignment'

    # create list of mismatch reads
    mismatch_list = df.loc[(df['match'] == 'mismatch') & (df['output_barcode'] != 'unclassified')]
    mismatch_list = mismatch_list[['filename','read_id', 'match', 'input_barcode', 'output_barcode']]

    # Save mismatch_list to a CSV file
    cleaned_parts = [part.replace("'", "") for part in file_name.split('_')[2:]]
    csv_name = f"mismatch_list_{'_'.join(cleaned_parts)}.csv"
    output_path = os.path.join(output_dir, csv_name)
    mismatch_list.to_csv(output_path, index=False)

    # Create a new DataFrame misclass_df with the relevant columns
    misclass_df = df[['input_barcode', 'output_barcode', 'match']]

    # Group by 'match' and calculate the count
    match_counts = misclass_df.groupby('match').size().reset_index(name=file_name.replace(".txt", ""))
    
    # Merge the current match counts with the combined match counts
    combined_match_counts = pd.merge(combined_match_counts, match_counts, 
                                     on='match', how='outer')

    return combined_match_counts

def main():
    parser = argparse.ArgumentParser(description='Process summary files and generate mismatch lists.')
    parser.add_argument('--input', required=True, help='Input directory path to dorado summary files')
    parser.add_argument('--output', required=True, help='Output directory path for 1) mismatch list')
    args = parser.parse_args()

    input_dir = args.input
    output_dir = args.output

    # Get the list of filenames in the input directory
    file_names = os.listdir(input_dir)

    # Initialize combined_match_counts
    combined_match_counts = pd.DataFrame(columns=['match'])

    for file_name in file_names:
        if file_name.endswith(".txt"):
            combined_match_counts = process_summary_file(input_dir, output_dir, file_name, combined_match_counts)

    # Save the combined match counts to a CSV file
    output_file = os.path.join(input_dir, 'match_counts_combined_1to1.csv')
    combined_match_counts.to_csv(output_file, index=False)

if __name__ == "__main__":
    main()
