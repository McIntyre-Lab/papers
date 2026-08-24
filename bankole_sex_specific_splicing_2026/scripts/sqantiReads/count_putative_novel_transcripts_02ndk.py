#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 11:38:31 2025

@author: nkeil
"""
import pandas as pd
import os

# Set paths to directories
rmg_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"

#Set paths to ujc counts file
mel_ujc_file = os.path.join(rmg_dir, "rmg_mel_comb_TRs_facet_sex_putative_underannotation.csv")
sim_ujc_file = os.path.join(rmg_dir, "rmg_sim_comb_TRs_facet_sex_putative_underannotation.csv")
yak_ujc_file = os.path.join(rlr_dir, "TO_dyak_comb_TRs_facet_sex_putative_underannotation.csv")
san_ujc_file = os.path.join(rlr_dir, "TO_dsan_comb_TRs_facet_sex_putative_underannotation.csv")
ser_ujc_file = os.path.join(axk_dir, "kopp_dser_comb_TRs_facet_sex_putative_underannotation.csv")

mel_gene_file = os.path.join(rmg_dir, "rmg_mel_comb_TRs_facet_sex_gene_classfication.csv")
sim_gene_file = os.path.join(rmg_dir, "rmg_sim_comb_TRs_facet_sex_gene_classfication.csv")
yak_gene_file = os.path.join(rlr_dir, "TO_dyak_comb_TRs_facet_sex_gene_classfication.csv")
san_gene_file = os.path.join(rlr_dir, "TO_dsan_comb_TRs_facet_sex_gene_classfication.csv")
ser_gene_file = os.path.join(axk_dir, "kopp_dser_comb_TRs_facet_sex_gene_classfication.csv")

#Make dictionary linking file
species_files = {
    "mel": {
        "ujc_file": mel_ujc_file,
        "gene_file": mel_gene_file,
        "tag": "dmel6"
    },
    "sim": {
        "ujc_file": sim_ujc_file,
        "gene_file": sim_gene_file,
        "tag": "dsim2"
    },
    "yak": {
        "ujc_file": yak_ujc_file,
        "gene_file": yak_gene_file,
        "tag": "dyak2"
    },
    "san": {
        "ujc_file": san_ujc_file,
        "gene_file": san_gene_file,
        "tag": "dsan1"
    },
    "ser": {
        "ujc_file": ser_ujc_file,
        "gene_file": ser_gene_file,
        "tag": "dser1"
    }
}

counts=[]
#categories =set()
transcripts_df = pd.DataFrame()


for species, files in species_files.items():
    ujc_file = files["ujc_file"]
    gene_file = files["gene_file"]
    tag = files["tag"]
    
    # Read the ujc file into a DataFrame
    ujc_df = pd.read_csv(ujc_file)
    
    #Subset to rows with putative novel transcript
    ujc_df = ujc_df[ujc_df['flag_putative_novel_transcript'] == 1]
    
    #Count number of genes and transcripts with putative novel transcript
    num_trancript = ujc_df['jxnHash'].nunique()
    num_NNC = len(ujc_df[ujc_df['structural_category'] == 'novel_not_in_catalog'])
    num_NIC = len(ujc_df[ujc_df['structural_category'] == 'novel_in_catalog'])
    num_gene = ujc_df['associated_gene'].nunique()
    
    # Store results in a dictionary
    counts.append({
        "species": species,
        "tag": tag,
        "num_novel_transcripts": num_trancript,
        "num_NIC": num_NIC,
        "num_NNC": num_NNC,
        "num_genes_w_novel_transcript": num_gene
    })
    
    ujc_df['species'] = species
    ujc_df['tag'] = tag
    transcripts_df = pd.concat([transcripts_df, ujc_df], ignore_index=True)
    

count_df = pd.DataFrame(counts)
count_df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/count_putative_novel_transcripts.csv', index = False)
    
transcripts_df.to_csv('/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/putative_novel_transcript_list.csv', index = False)    