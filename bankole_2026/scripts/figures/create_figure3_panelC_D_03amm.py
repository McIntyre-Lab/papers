#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
same as B except restricted to genes 
with ge 2 full length transcripts (nonISM)

multiexon and multi transcript
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY = f"{PROJ}/zenodo/summary_files"
OUT = f"{PROJ}/Figures_amm"
os.makedirs(OUT, exist_ok=True)

species_list = ["dmel", "dsim", "dyak", "dsan", "dser"]

dark = {"dmel": "#966729", "dsim": "#3F78C1", "dyak": "#717273", "dsan": "#28827A", "dser": "#825CA6"}
medium = {"dmel": "#D4B08C", "dsim": "#9FC1E8", "dyak": "#B8B9BA", "dsan": "#8FCCC7", "dser": "#C3A9D3"}
light = {"dmel": "#F4D8A8", "dsim": "#C8E2F8", "dyak": "#D9DADB", "dsan": "#C6E9E5", "dser": "#E1D8F0"}

panel_c_cols = ["alt_IR", "alt_donor_acceptor"]
panel_c_labels = ["Alt. intron retention", "Alt. donor/acceptor"]

panel_d_cols = ["alt_5p_ER", "alt_3p_ER", "alt_ERSkip"]
panel_d_labels = ["Alt. 5′ ER", "Alt. 3′ ER", "Alt. ER skip"]

all_alt_cols = panel_c_cols + panel_d_cols
plot_data = {col: {} for col in all_alt_cols}
gene_counts = {species: {"panel_c": 0, "panel_d": 0} for species in species_list}

for species in species_list:
    gene_df = pd.read_csv(
        f"{SUMMARY}/gene_summary_min10reads_{species}.csv",
        usecols=["geneID", "numExon_GM", "num_ujc", "num_ism_ujc"],
        low_memory=False
    )

    as_df = pd.read_csv(
        f"{SUMMARY}/gene_summary_min10reads_from_AS_analysis_{species}.csv",
        usecols=["geneID"] + all_alt_cols,
        low_memory=False
    ).drop_duplicates("geneID")

    df = gene_df.merge(as_df, on="geneID", how="left", validate="one_to_one")

    for col in all_alt_cols:
        df[col] = df[col].fillna("none")

    df["num_nonISM_UJC"] = df["num_ujc"] - df["num_ism_ujc"]

    df = df[
        (df["numExon_GM"] >= 2)
        & (df["num_nonISM_UJC"] >= 2)
    ].copy()

    ## count genes going into C and D
    panel_c_mask = (df[panel_c_cols] != "none").any(axis=1)
    gene_counts[species]["panel_c"] = panel_c_mask.sum()

    panel_d_mask = (df[panel_d_cols] != "none").any(axis=1)
    gene_counts[species]["panel_d"] = panel_d_mask.sum()

    
    denominator = len(df)
    print(f"total genes for {species} is {denominator}")
    
    for col in all_alt_cols:
        if denominator == 0:
            plot_data[col][species] = {"anno_UJC": 0, "anno_ERP": 0, "novel_ERP": 0}
        else:
            plot_data[col][species] = {
                "anno_UJC": (df[col] == "anno_UJC").sum() / denominator,
                "anno_ERP": (df[col] == "anno_ERP").sum() / denominator,
                "novel_ERP": (df[col] == "novel_ERP").sum() / denominator,
            }

for species in species_list:
    print(f"{species}:  panel C = {gene_counts[species]['panel_c']},  panel D = {gene_counts[species]['panel_d']}")


def make_panel(cols, labels, outfile):
    fig, ax = plt.subplots(figsize=(11, 7))
    x = np.arange(len(cols))
    bar_width = 0.14
    offsets = np.arange(len(species_list)) - (len(species_list) - 1) / 2

    for i, species in enumerate(species_list):
        x_pos = x + offsets[i] * bar_width
        anno_ujc = np.array([plot_data[col][species]["anno_UJC"] for col in cols])
        anno_erp = np.array([plot_data[col][species]["anno_ERP"] for col in cols])
        novel_erp = np.array([plot_data[col][species]["novel_ERP"] for col in cols])

        ax.bar(x_pos, anno_ujc, bar_width, color=dark[species], edgecolor="black", linewidth=0.5)
        ax.bar(x_pos, anno_erp, bar_width, bottom=anno_ujc, color=medium[species], edgecolor="black", linewidth=0.5)
        ax.bar(x_pos, novel_erp, bar_width, bottom=anno_ujc + anno_erp, color=light[species], edgecolor="black", linewidth=0.5)

    ax.set_ylabel("Proportion of genes", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    plt.tight_layout()
    plt.savefig(outfile, format="pdf", bbox_inches="tight", dpi=600)
    plt.close()

make_panel(panel_c_cols, panel_c_labels, f"{OUT}/figure3_panelC_02mdg.pdf")
make_panel(panel_d_cols, panel_d_labels, f"{OUT}/figure3_panelD_02mdg.pdf")
