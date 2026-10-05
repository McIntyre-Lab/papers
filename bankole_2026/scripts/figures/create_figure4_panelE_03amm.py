#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
full length (no ISM)
multiexon AND multitranscript

"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
SUMMARY = f"{PROJ}/zenodo/summary_files"
OUT = f"{PROJ}/Figures_amm"
os.makedirs(OUT, exist_ok=True)

species_list = ["dmel", "dsim", "dser"]

dark = {"dmel": "#966729", "dsim": "#3F78C1", "dyak": "#717273", "dsan": "#28827A", "dser": "#825CA6"}
medium = {"dmel": "#D4B08C", "dsim": "#9FC1E8", "dyak": "#B8B9BA", "dsan": "#8FCCC7", "dser": "#C3A9D3"}
light = {"dmel": "#F4D8A8", "dsim": "#C8E2F8", "dyak": "#D9DADB", "dsan": "#C6E9E5", "dser": "#E1D8F0"}

same_erp_cols = ["alt_IR", "alt_donor_acceptor"]
same_erp_labels = ["Alt. intron retention", "Alt. donor/acceptor"]

different_erp_cols = ["alt_5p_ER", "alt_3p_ER", "alt_ERSkip"]
different_erp_labels = ["Alt. 5′ ER", "Alt. 3′ ER", "Alt. ER skip"]

all_alt_cols = same_erp_cols + different_erp_cols
plot_data = {col: {} for col in all_alt_cols}

for species in species_list:
    gene_df = pd.read_csv(
        f"{SUMMARY}/gene_summary_min10reads_{species}.csv",
        usecols=[
            "geneID",
            "numExon_GM",
            "num_ujc",
            "num_ism_ujc",
            "num_ujc_bias_F",
            "num_ujc_bias_M",
            "num_ERPp_bias_F",
            "num_ERPp_bias_M",
        ],
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

    ## define denominator - multiexon genes with multiple transcripts that show sex bias
    sex_biased = (
        df["num_ujc_bias_F"] + df["num_ujc_bias_M"]
        + df["num_ERPp_bias_F"] + df["num_ERPp_bias_M"]
    ) >= 1

    df = df[
        (df["numExon_GM"] >= 2)
        & (df["num_nonISM_UJC"] >= 2)
        & sex_biased
    ].copy()

    denominator = len(df)
    print(f'genes with sex bias and multiple observed transcripts with at least 10 reads in {species}:')
    print(denominator)

    # calc proportions
    for col in all_alt_cols:
        plot_data[col][species] = {
            "anno_UJC": (df[col] == "anno_UJC").sum() / denominator,
            "anno_ERP": (df[col] == "anno_ERP").sum() / denominator,
            "novel_ERP": (df[col] == "novel_ERP").sum() / denominator,
        }

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

    ax.set_ylabel("Proportion of sex-biased genes", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    plt.tight_layout()
    plt.savefig(outfile, format="pdf", bbox_inches="tight", dpi=600)
    plt.close()

make_panel(same_erp_cols, same_erp_labels, f"{OUT}/figure4_panelE_sameERP_02mdg.pdf")
make_panel(different_erp_cols, different_erp_labels, f"{OUT}/figure4_panelE_differentERP_02mdg.pdf")
