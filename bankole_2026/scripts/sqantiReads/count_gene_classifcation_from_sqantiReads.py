#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 23:17:46 2025

@author: nkeil
"""
import os
import pandas as pd

# Directory with gene summary files
ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables"

# Hard-coded classification rules
classification_rules = [
    {
        "class_name": "annotated_with_well_covered_FSM",
        "flag_ge_100_reads": 1,
        "flag_has_FSM_ujc": 1,
        "flag_has_FSM_ujc_ge_20_perc": 1
    },
    {
        "class_name": "annotated_with_low_coverage_FSM",
        "flag_ge_100_reads": 1,
        "flag_has_FSM_ujc": 1,
        "flag_has_FSM_ujc_ge_20_perc": 0
    },
    {
        "class_name": "underannotated_no_candidate_transcripts",
        "flag_ge_100_reads": 1,
        "flag_has_FSM_ujc": 0,
        "flag_has_any_ujc_ge_20_perc": 0
    },
    {
        "class_name": "underannotated_with_candidate_transcript",
        "flag_ge_100_reads": 1,
        "flag_has_FSM_ujc": 0,
        "flag_has_any_ujc_ge_20_perc": 1
    },
    {
        "class_name": "ge_85_perc_w_anno_transcript",
        "flag_ge_100_reads": 1,
        "flag_gt_85perc_w_anno_transcript": 1
    }
]

# Store results
all_counts = []

for file in os.listdir(ind):
    if not file.endswith("_gene_ujc_summary_from_sqantiReads.csv"):
        continue

    species = file.replace("_gene_ujc_summary_from_sqantiReads.csv", "")
    gene_df = pd.read_csv(os.path.join(ind, file))
    
    gene_df["flag_ge_100_reads"] = (gene_df["total_reads"] >= 100).astype(int)


    # Count thresholds
    counts = {
        "species": species,
        "ge_1_reads": (gene_df["total_reads"] >= 1).sum(),
        "ge_25_reads": (gene_df["total_reads"] >= 25).sum(),
        "ge_50_reads": (gene_df["total_reads"] >= 50).sum(),
        "ge_100_reads": (gene_df["total_reads"] >= 100).sum()
    }
    
    

    # Apply each rule
    for rule in classification_rules:
        rule_name = rule["class_name"]
        mask = pd.Series([True] * len(gene_df))
        for flag, value in rule.items():
            if flag == "class_name":
                continue
            mask &= (gene_df[flag] == value)
        counts[rule_name] = mask.sum()

    all_counts.append(counts)

# Final table
summary_df = pd.DataFrame(all_counts)
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
summary_df.to_csv(os.path.join(outd,"count_gene_classification.csv"), index=False)
