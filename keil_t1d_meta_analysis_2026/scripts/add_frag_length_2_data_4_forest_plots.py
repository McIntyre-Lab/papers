#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 29 13:44:52 2025

@author: nkeil
"""
import os
import pandas as pd

# Base path and list of subdirectories
base_path = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts"
subdirs = ["data_4_forest_plots_ES_diff", "data_4_forest_plots", "data_4_forest_plots_means"]  # replace with actual subdirectory names

for subdir in subdirs:
    dir_path = os.path.join(base_path, subdir)

    # Skip if the directory doesn't exist
    if not os.path.isdir(dir_path):
        print(f"Directory not found: {dir_path}")
        continue

    # Loop over files in the subdirectory
    for filename in os.listdir(dir_path):
        if filename.endswith("features_FE.csv"):
            file_path = os.path.join(dir_path, filename)
            print(f"Processing {file_path}...")

            # Read and modify the CSV
            df = pd.read_csv(file_path)
            df["frag_length"] = df["ef_end"] - df["ef_start"]

            # Save to a new file in the same directory
            output_file = filename.replace("features_FE.csv", "features_FE_with_fraglen.csv")
            df.to_csv(os.path.join(dir_path, output_file), index=False)

            print(f"Saved to {output_file}")