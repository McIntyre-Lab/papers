#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 21:11:56 2025

@author: nkeil
"""
import pandas as pd
import os

# Set paths to directories
rmg_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"

# Set paths to gene counts files
mel_gene_file = os.path.join(rmg_dir, "rmg_mel_comb_TRs_facet_sex_gene_counts.csv")
sim_gene_file = os.path.join(rmg_dir, "rmg_sim_comb_TRs_facet_sex_gene_counts.csv")
yak_gene_file = os.path.join(rlr_dir, "TO_dyak_comb_TRs_facet_sex_gene_counts.csv")
san_gene_file = os.path.join(rlr_dir, "TO_dsan_comb_TRs_facet_sex_gene_counts.csv")
ser_gene_file = os.path.join(axk_dir, "kopp_dser_comb_TRs_facet_sex_gene_counts.csv")

#Set paths to ujc counts file
mel_ujc_file = os.path.join(rmg_dir, "rmg_mel_comb_TRs_facet_sex_ujc_counts.csv")
sim_ujc_file = os.path.join(rmg_dir, "rmg_sim_comb_TRs_facet_sex_ujc_counts.csv")
yak_ujc_file = os.path.join(rlr_dir, "TO_dyak_comb_TRs_facet_sex_ujc_counts.csv")
san_ujc_file = os.path.join(rlr_dir, "TO_dsan_comb_TRs_facet_sex_ujc_counts.csv")
ser_ujc_file = os.path.join(axk_dir, "kopp_dser_comb_TRs_facet_sex_ujc_counts.csv")

#Make dictionary linking file
species_files = {
    "mel": {
        "gene_file": mel_gene_file,
        "ujc_file": mel_ujc_file,
        "tag": "dmel6"
    },
    "sim": {
        "gene_file": sim_gene_file,
        "ujc_file": sim_ujc_file,
        "tag": "dsim2"
    },
    "yak": {
        "gene_file": yak_gene_file,
        "ujc_file": yak_ujc_file,
        "tag": "dyak2"
    },
    "san": {
        "gene_file": san_gene_file,
        "ujc_file": san_ujc_file,
        "tag": "dsan1"
    },
    "ser": {
        "gene_file": ser_gene_file,
        "ujc_file": ser_ujc_file,
        "tag": "dser1"  # No tag provided for "ser"
    }
}

for species, files in species_files.items():
    gene_file = files["gene_file"]
    tag = files["tag"]
    
    # Read the gene file into a DataFrame
    gene_df = pd.read_csv(gene_file)
    
    # Filter for annotated genes
    anno_gene_df = gene_df[gene_df["flag_annotated_gene"] == 1]
    
    # Group by associated_gene and sum the total_read_count across all samples
    gene_sum_df = anno_gene_df.groupby("associated_gene", as_index=False)[
    ["total_read_count", "full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "unique_jxnHash_counts"]].sum()

    gene_sum_df = gene_sum_df.rename(columns={
    "full-splice_match": "FSM_reads",
    "incomplete-splice_match": "ISM_reads",
    "novel_in_catalog": "NIC_reads",
    "novel_not_in_catalog": "NNC_reads",
    "total_read_count": "total_reads",
    "num_jxnHash": "unique_jxnHash_counts"
        })

    #  Add percentage columns
    for col in ["FSM_reads", "ISM_reads", "NIC_reads", "NNC_reads"]:
        abbr = col.split("_")[0]  # e.g., FSM
        gene_sum_df[f"perc_reads_{abbr}"] = (
            gene_sum_df[col] / gene_sum_df["total_reads"] * 100
            ).round(5)  # Rounded for neatness
    
    #Sum perc FSM+perc ISM + perc NIC
    gene_sum_df["perc_w_anno_transcript"] = gene_sum_df["perc_reads_FSM"] + gene_sum_df["perc_reads_ISM"] + gene_sum_df["perc_reads_NIC"]
    
    #Get ujc info
    ujc_file = files["ujc_file"]
    ujc_df = pd.read_csv(ujc_file)

     # Filter for annotated genes
    anno_ujc_df = ujc_df[ujc_df["flag_annotated_gene"] == 1]
    
    #Structural categories of interest
    # Define only the categories you care about
    structural_categories = [
        "full-splice_match",
        "incomplete-splice_match",
        "novel_in_catalog",
        "novel_not_in_catalog"
    ]

    # Subset to those rows only
    ujc_filt_df = anno_ujc_df[anno_ujc_df["structural_category"].isin(structural_categories)].copy()
    
    # Map long-form to abbreviation
    abbr_map = {
        "full-splice_match": "FSM",
        "incomplete-splice_match": "ISM",
        "novel_in_catalog": "NIC",
        "novel_not_in_catalog": "NNC"
    }
    
    ujc_filt_df["structural_abbr"] = ujc_filt_df["structural_category"].map(abbr_map)

    # Group and pivot
    ujc_cnt_df = ujc_filt_df.groupby(["associated_gene", "structural_abbr"])["jxnHash"].nunique().reset_index()
    ujc_cnt_df = ujc_cnt_df.pivot(
        index="associated_gene",
        columns="structural_abbr",
        values="jxnHash"
    ).fillna(0).astype(int)
    
    # Rename columns
    ujc_cnt_df.columns = [f"{col}_ujcs" for col in ujc_cnt_df.columns]
    ujc_cnt_df = ujc_cnt_df.reset_index()
    
    #Merge ujc counts and NIC counts
    gene_ujc_cnt_df = pd.merge(gene_sum_df, ujc_cnt_df, on="associated_gene", how="outer", indicator=True)
    
    
    mask = gene_ujc_cnt_df["_merge"] == "left_only"
    numeric_cols = gene_ujc_cnt_df.select_dtypes(include=["number"]).columns
    gene_ujc_cnt_df.loc[mask, numeric_cols] = gene_ujc_cnt_df.loc[mask, numeric_cols].fillna(0)
    
    gene_ujc_cnt_df = gene_ujc_cnt_df.drop(columns=["_merge"])
    
    #Mkae flags
    gene_ujc_cnt_df["flag_ge_100_reads"] = (gene_ujc_cnt_df["total_reads"] >= 100).astype(int)
    gene_ujc_cnt_df["flag_has_FSM_ujc"] = (gene_ujc_cnt_df["FSM_ujcs"] > 0).astype(int)
    gene_ujc_cnt_df["flag_gt_85perc_w_anno_transcript"] = (gene_ujc_cnt_df["perc_w_anno_transcript"] > 85).astype(int)

    
   # Sum read counts per gene per jxnHash
    ujc_sum_df = ujc_filt_df.groupby(["associated_gene", "jxnHash", "structural_category"])["read_count"].sum().reset_index()

    #  Merge with gene-level total_reads
    ujc_sum_df = ujc_sum_df.merge(
        gene_sum_df[["associated_gene", "total_reads"]],
        on="associated_gene",
        how="left"
    )

    #Calculate percent of gene reads per jxnHash
    ujc_sum_df["perc_of_gene_reads"] = ujc_sum_df["read_count"] / ujc_sum_df["total_reads"] * 100

    # Create flags at jxnHash level
    ujc_sum_df["is_FSM_ge_20"] = (ujc_sum_df["structural_category"] == "full-splice_match") & (ujc_sum_df["perc_of_gene_reads"] >= 20)
    ujc_sum_df["is_any_ujc_ge_20"] = ujc_sum_df["perc_of_gene_reads"] >= 20

    # Aggregate to gene-level flags
    fsm_flag_df = ujc_sum_df.groupby("associated_gene")["is_FSM_ge_20"].any().reset_index()
    fsm_flag_df = fsm_flag_df.rename(columns={"is_FSM_ge_20": "flag_has_FSM_ujc_ge_20_perc"})

    any_flag_df = ujc_sum_df.groupby("associated_gene")["is_any_ujc_ge_20"].any().reset_index()
    any_flag_df = any_flag_df.rename(columns={"is_any_ujc_ge_20": "flag_has_any_ujc_ge_20_perc"})

    # Merge flags into main gene-level dataframe
    gene_ujc_cnt_df =   gene_ujc_cnt_df.merge(fsm_flag_df, on="associated_gene", how="left")
    gene_ujc_cnt_df =   gene_ujc_cnt_df.merge(any_flag_df, on="associated_gene", how="left")
    
    # Convert booleans to 0/1
    gene_ujc_cnt_df["flag_has_FSM_ujc_ge_20_perc"] = gene_ujc_cnt_df["flag_has_FSM_ujc_ge_20_perc"].fillna(False).astype(int)
    gene_ujc_cnt_df["flag_has_any_ujc_ge_20_perc"] = gene_ujc_cnt_df["flag_has_any_ujc_ge_20_perc"].fillna(False).astype(int)


    out = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables"
    output_file = os.path.join(out, f"{tag}_gene_ujc_summary_from_sqantiReads.csv")
    gene_ujc_cnt_df.to_csv(output_file, index=False)
    
    
    