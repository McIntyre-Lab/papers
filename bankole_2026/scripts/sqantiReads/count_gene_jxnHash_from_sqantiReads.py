#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 14:06:28 2025

@author: nkeil
"""
import pandas as pd
import os

# Set paths to directories
rmg_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"

anno_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary"


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

#Get gene counts for all 5 species
gene_counts = []

for species, files in species_files.items():
    gene_file = files["gene_file"]
    tag = files["tag"]
    
    # Read the gene file into a DataFrame
    gene_df = pd.read_csv(gene_file)
    
    # Filter for annotated genes
    anno_gene_df = gene_df[gene_df["flag_annotated_gene"] == 1]
    
    # Group by associated_gene and sum the total_read_count across all samples
    gene_sum_df = anno_gene_df.groupby("associated_gene", as_index=False)["total_read_count"].sum()
    
    # Now apply the thresholds
    num_gene_ge1_read = (gene_sum_df["total_read_count"] >= 1).sum()
    num_gene_ge25_read = (gene_sum_df["total_read_count"] >= 25).sum()
    num_gene_ge50_read = (gene_sum_df["total_read_count"] >= 50).sum()
    num_gene_ge100_read = (gene_sum_df["total_read_count"] >= 100).sum()
    
    # Store results in a dictionary
    gene_counts.append({
        "species": species,
        "tag": tag,
        "num_gene_ge1_read": num_gene_ge1_read,
        "num_gene_ge25_read": num_gene_ge25_read,
        "num_gene_ge50_read": num_gene_ge50_read,
        "num_gene_ge100_read": num_gene_ge100_read
    })

# Create a df of counts
gene_count_df = pd.DataFrame(gene_counts)

#Get ujc_counts for all 5 species
ujc_counts = []
anno_ujc_counts = []

for species, files in species_files.items():
    ujc_file = files["ujc_file"]
    tag = files["tag"]
    
    # Read the ujc file into a DataFrame
    ujc_df = pd.read_csv(ujc_file)    
    
    # Group by jxnHash and sum the total_read_count across all samples
    ujc_sum_df = ujc_df.groupby("jxnHash", as_index=False)["read_count"].sum()
    
    # Now apply the thresholds
    num_ujc_ge1_read = (ujc_sum_df["read_count"] >= 1).sum()
    num_ujc_ge25_read = (ujc_sum_df["read_count"] >= 25).sum()
    num_ujc_ge50_read = (ujc_sum_df["read_count"] >= 50).sum()
    num_ujc_ge100_read = (ujc_sum_df["read_count"] >= 100).sum()
    
    # Store results in a dictionary
    ujc_counts.append({
        "species": species,
        "tag": tag,
        "num_jxnHash_ge1_read": num_ujc_ge1_read,
        "num_jxnHash_ge25_read": num_ujc_ge25_read,
        "num_jxnHash_ge50_read": num_ujc_ge50_read,
        "num_jxnHash_ge100_read": num_ujc_ge100_read
    })
    
# =============================================================================
#     anno_file_name = f"fiveSpecies_{tag}_full_annotation_w_dataFlag.csv"
#     anno_col = f"{tag}_jxnHash"
#     anno_file = os.path.join(anno_dir, anno_file_name)
#     
#     anno_df = pd.read_csv(anno_file,usecols=[anno_col] )
#     anno_ujc_df = ujc_df[ujc_df["jxnHash"].isin(anno_df[anno_col])]
# =============================================================================
    
ujc_count_df = pd.DataFrame(ujc_counts)

#merge with indicator check
df = pd.merge(gene_count_df, ujc_count_df, on=['species','tag'], how='outer', indicator=True)

df = df[df['_merge'] == 'both'].drop(columns=['_merge'])

df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/gene_ujc_counts_from_sqantiReads.csv')

