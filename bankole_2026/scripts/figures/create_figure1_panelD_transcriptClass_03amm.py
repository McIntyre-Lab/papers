#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Set directories, species order, and plotting colors
PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY = f"{PROJ}/zenodo/summary_files"
OUT = f"{PROJ}/Figures_amm"

species_list = ["dmel", "dsim", "dyak", "dsan", "dser"]

species_colors_detected = {"dmel": "#966729", "dsim": "#3F78C1", "dyak": "#717273", "dsan": "#28827A", "dser": "#825CA6"}
species_colors_undetected = {"dmel": "#D4B08C", "dsim": "#9FC1E8", "dyak": "#B8B9BA", "dsan": "#8FCCC7", "dser": "#C3A9D3"}

# Read component summary and prepare gene-level information for each species
comp_df = pd.read_csv(f"{SUMMARY}/component_summary.csv", low_memory=False)
species_results = {}

for species in species_list:
    gene_col = f"geneID_concat_{species}"
    simple_flag = f"flag_simple_{species}"
    monoExon_flag = f"flag_monoExon_{species}"

    # Expand concatenated gene IDs into gene/component rows
    gene_comp = comp_df[["componentID", gene_col, simple_flag, monoExon_flag]].dropna(subset=[gene_col]).copy()
    gene_comp["geneID"] = gene_comp[gene_col].str.split("|")
    gene_comp = gene_comp.explode("geneID")
    gene_comp["geneID"] = gene_comp["geneID"].str.strip()

    # Count components and summarize simple-component and monoexon flags per gene
    gene_df = gene_comp.groupby("geneID").agg(
        num_components=("componentID", "nunique"),
        num_simple_components=(simple_flag, "sum"),
        flag_monoExon=(monoExon_flag, "min"),
    ).reset_index()
    
    ## flag_monoExon = 1 only if all components are monoexon?  yes
    # mixed = gene_comp.groupby("geneID")[monoExon_flag].nunique()
    # mixed = mixed[mixed == 2].index
    # gene_mixed = gene_comp[
    #     gene_comp["geneID"].isin(mixed)].sort_values("geneID")    
    # find = gene_df[gene_df["geneID"] == "GFF5_"]
    
    # mono_genes = gene_comp.groupby("geneID")[monoExon_flag].min()
    # mono_genes = mono_genes[mono_genes == 1].index

    # gene_comp_mono = gene_comp[
    # gene_comp["geneID"].isin(mono_genes)
    # ]



    # Assign transcript class, with monoexon taking precedence over single transcript
    gene_df["transcriptClass"] = "multiTranscript"
    single_condition = (gene_df["num_components"] == 1) & (gene_df["num_simple_components"] == 1)
    gene_df.loc[single_condition, "transcriptClass"] = "singleTranscript"
    gene_df.loc[gene_df["flag_monoExon"] == 1, "transcriptClass"] = "monoExon"

    # Flag genes present in the gene summary as detected and identify GFF genes
    summary_df = pd.read_csv(f"{SUMMARY}/gene_summary_{species}.csv", usecols=["geneID"], low_memory=False)
    gene_df["flag_detected"] = gene_df["geneID"].isin(summary_df["geneID"]).astype(int)
    gene_df["flag_GFF"] = gene_df["geneID"].str.startswith("GFF").astype(int)
    species_results[species] = gene_df

# Define the two plots and transcript categories
plot_groups = {"GFF": 1, "nonGFF": 0}
categories = ["monoExon", "singleTranscript", "multiTranscript"]
x = np.arange(len(categories))
bar_width = 0.15
offsets = np.arange(len(species_list)) - (len(species_list) - 1) / 2

# Loop over gene groups and species
for group_name, flag_gff in plot_groups.items():
    fig, ax = plt.subplots(figsize=(10, 7))
    for i, species in enumerate(species_list):
        df = species_results[species]
        conditions = {
            "monoExon": (df["flag_GFF"] == flag_gff) & (df["transcriptClass"] == "monoExon"),
            "singleTranscript": (df["flag_GFF"] == flag_gff) & (df["transcriptClass"] == "singleTranscript"),
            "multiTranscript": (df["flag_GFF"] == flag_gff) & (df["transcriptClass"] == "multiTranscript"),
        }

 
        # Count detected and undetected genes for each condition
        detected, undetected = [], []
        for category, condition in conditions.items():
            detected.append((condition & (df["flag_detected"] == 1)).sum())
            undetected.append((condition & (df["flag_detected"] == 0)).sum())

        ## amm added counts below 
        totals = np.array(detected) + np.array(undetected)
        overall_total = int(np.sum(totals))
        for category, total in zip(categories, totals):
            print(group_name, species, category, int(total))
        print(f"{group_name} | {species} | overall total: {overall_total}")
            
        # Stack undetected genes above detected genes using species-specific colors
        xpos = x + offsets[i] * bar_width
        ax.bar(xpos, detected, bar_width, color=species_colors_detected[species],
               edgecolor="black", linewidth=0.5, label=species)
        ax.bar(xpos, undetected, bar_width, bottom=detected, color=species_colors_undetected[species],
               edgecolor="black", linewidth=0.5)
        
        # amm set y axis
        ax.set_ylim(0, 12000)
        
    # Label axes and save one PDF for each gene group
    ax.set_xlabel("Transcript Class")
    ax.set_ylabel("Number of Genes")
    ax.set_xticks(x)
    ax.set_xticklabels(["Mono-exon", "Single transcript", "Multi-transcript"])
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/figure1_panelD_{group_name}_transcriptClass_02mdg.pdf", bbox_inches="tight")
    plt.close(fig)