#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Set directories, species order, and plotting colors
PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO = f"{PROJ}/zenodo"
SUMMARY = f"{ZENODO}/summary_files"
OUT = f"{PROJ}/Figures_amm"

species_list = ["dmel", "dsim", "dyak", "dsan", "dser"]
species_2_genome = {"dmel": "dmel6", "dsim": "dsim2", "dyak": "dyak2", "dsan": "dsan1", "dser": "dser1"}

species_colors_detected = {"dmel": "#966729", "dsim": "#3F78C1", "dyak": "#717273", "dsan": "#28827A", "dser": "#825CA6"}
species_colors_undetected = {"dmel": "#D4B08C", "dsim": "#9FC1E8", "dyak": "#B8B9BA", "dsan": "#8FCCC7", "dser": "#C3A9D3"}

# Read annotations and flag detected genes, gene-model exon-region counts, and GFF genes
species_results = {}

for species in species_list:
    anno_df = pd.read_csv(f"{ZENODO}/fiveSpecies_{species_2_genome[species]}_full_annotation.csv", low_memory=False)
    summary_df = pd.read_csv(f"{SUMMARY}/gene_summary_{species}.csv", usecols=["geneID"], low_memory=False)

    gene_df = anno_df.groupby("geneID")["ERP"].first().reset_index()
    gene_df["flag_detected"] = gene_df["geneID"].isin(summary_df["geneID"]).astype(int)
    gene_df["numExon_GM"] = gene_df["ERP"].str.len() - 2
    gene_df["flag_GFF"] = gene_df["geneID"].str.startswith("GFF").astype(int)
    species_results[species] = gene_df

# Define the two plots and exon-region categories
plot_groups = {"GFF": 1, "nonGFF": 0}
categories = ["1", "2", "3", "4", "5+"]
x = np.arange(len(categories))
bar_width = 0.15
offsets = np.arange(len(species_list)) - (len(species_list) - 1) / 2

# Loop over gene groups and species
for group_name, flag_gff in plot_groups.items():
    fig, ax = plt.subplots(figsize=(10, 7))
    for i, species in enumerate(species_list):
        df = species_results[species]
        conditions = {
            "1": (df["flag_GFF"] == flag_gff) & (df["numExon_GM"] == 1),
            "2": (df["flag_GFF"] == flag_gff) & (df["numExon_GM"] == 2),
            "3": (df["flag_GFF"] == flag_gff) & (df["numExon_GM"] == 3),
            "4": (df["flag_GFF"] == flag_gff) & (df["numExon_GM"] == 4),
            "5+": (df["flag_GFF"] == flag_gff) & (df["numExon_GM"] >= 5),
        }

        # Count detected and undetected genes for each condition
        detected, undetected = [], []
        for category, condition in conditions.items():
            det = (condition & (df["flag_detected"] == 1)).sum()
            undet = (condition & (df["flag_detected"] == 0)).sum()
            detected.append(det)
            undetected.append(undet)
        print(f"{group_name} | {species} | total genes = {sum(detected) + sum(undetected)}")
            
                    
        # Stack undetected genes above detected genes using species-specific colors
        xpos = x + offsets[i] * bar_width
        ax.bar(xpos, detected, bar_width, color=species_colors_detected[species],
               edgecolor="black", linewidth=0.5, label=species)
        ax.bar(xpos, undetected, bar_width, bottom=detected, color=species_colors_undetected[species],
               edgecolor="black", linewidth=0.5)
        ## amm added setting yaxis
        ax.set_ylim(0, 6000)
        
    # Label axes and save one PDF for each gene group
    ax.set_xlabel("Number of Exon Regions in Gene Model")
    ax.set_ylabel("Number of Genes")
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{OUT}/figure1_panelD_{group_name}_numExon_02mdg.pdf", bbox_inches="tight")
    plt.close(fig)