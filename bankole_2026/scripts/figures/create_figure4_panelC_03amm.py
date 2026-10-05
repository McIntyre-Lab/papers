#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
among orthologous genes that can be analyzed in all species,
 what fraction show the same sex-bias direction conserved across species?

"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BASE_DIR = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures_amm"
SUMMARY_DIR = f"{BASE_DIR}/zenodo/summary_files"
TABLES_DIR = f"{BASE_DIR}/Tables"
MERGED_ORTHOLOG_FILE = f"{TABLES_DIR}/gene_sexBias_merged_on_genesetid.csv"

sex_bias_species_list = ["dmel", "dsim", "dser"]
analyzed_categories = ["F", "M", "B", "unbiased"]

sexbias_color = {
    "F": "#FF0000",
    "M": "#0000FF",
    "B": "#A020F0",
}


# Per-species sex-bias proportions
sexbias_plot_data = {}

for species in sex_bias_species_list:
    df = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{species}.csv", low_memory=False)

    # create True for each row where sexbias F,M,B or unbiased, else F
    analyzable = df["sexBias"].isin(analyzed_categories)
    n = analyzable.sum()  # count T

    print(f"\n{species}")
    print("denominator:", n)
    print("M:", (analyzable & (df["sexBias"] == "M")).sum())
    print("F:", (analyzable & (df["sexBias"] == "F")).sum())
    print("B:", (analyzable & (df["sexBias"] == "B")).sum())
    print("U:", (analyzable & (df["sexBias"] == "unbiased")).sum())

    # M, F, and B proportions
    sexbias_plot_data[species] = {
        "M_prop": (analyzable & (df["sexBias"] == "M")).sum() / n if n > 0 else 0,
        "F_prop": (analyzable & (df["sexBias"] == "F")).sum() / n if n > 0 else 0,
        "B_prop": (analyzable & (df["sexBias"] == "B")).sum() / n if n > 0 else 0,
    }


# Conserved sex bias across genesets
merged_df = pd.read_csv(MERGED_ORTHOLOG_FILE, low_memory=False)

comparisons = [
    ("dmel and dser", ["dmel", "dser"]),
    ("dmel and dsim", ["dmel", "dsim"]),
    ("dser and dsim", ["dser", "dsim"]),
    ("dmel, dsim, and dser", ["dmel", "dsim", "dser"]),
]

gap = 0.3
bars = [(i + 1, sexbias_plot_data[sp]) for i, sp in enumerate(sex_bias_species_list)]
comp_x0 = len(sex_bias_species_list) + 1 + gap
comp_labels = []

for j, (label, species_group) in enumerate(comparisons):

    # Genesets analyzable in every species
    # keep orthologs with at least 1 analyzable gene in every species in group
    denom = pd.Series(True, index=merged_df.index)
    for sp in species_group:
        denom &= merged_df[f"{sp}_num_analyzable"] >= 1

    n = denom.sum()  # count orthologs passing filter

    # Conserved M, F, and B bias
    # ortholog conserved F/M/B bias only if every species in group has at least 1 F/M/B biased gene
    conserved_M = pd.Series(True, index=merged_df.index)
    conserved_F = pd.Series(True, index=merged_df.index)
    conserved_B = pd.Series(True, index=merged_df.index)

    for sp in species_group:
        conserved_M &= merged_df[f"{sp}_num_M"] >= 1
        conserved_F &= merged_df[f"{sp}_num_F"] >= 1
        conserved_B &= merged_df[f"{sp}_num_B"] >= 1

    print(f"\n{label}")
    print("denominator (orthologs):", n)
    print("conserved M:", (denom & conserved_M).sum())
    print("conserved F:", (denom & conserved_F).sum())
    print("conserved B:", (denom & conserved_B).sum())
    


    # count orthologs analyzable across all and conserved in bias, divide by n
    props = {
        "M_prop": (denom & conserved_M).sum() / n if n > 0 else 0,
        "F_prop": (denom & conserved_F).sum() / n if n > 0 else 0,
        "B_prop": (denom & conserved_B).sum() / n if n > 0 else 0,
    }

    bars.append((comp_x0 + j, props))
    comp_labels.append(label)


# Plot stacked proportions
fig, ax = plt.subplots(figsize=(16, 7))
bar_width = 0.5

for xpos, props in bars:
    ax.bar(xpos, props["M_prop"], bar_width, color=sexbias_color["M"], edgecolor="black", linewidth=0.5)
    ax.bar(xpos, props["F_prop"], bar_width, bottom=props["M_prop"], color=sexbias_color["F"], edgecolor="black", linewidth=0.5)
    ax.bar(xpos, props["B_prop"], bar_width, bottom=props["M_prop"] + props["F_prop"], color=sexbias_color["B"], edgecolor="black", linewidth=0.5)

# Separate species and conserved comparisons
ax.axvline(x=(len(sex_bias_species_list) + comp_x0) / 2, color="gray", linestyle="--", linewidth=1, alpha=0.6)

# Axis labels
ax.set_xticks([xpos for xpos, _ in bars])
ax.set_xticklabels(sex_bias_species_list + comp_labels, fontsize=9)
ax.set_xlabel("Species", fontsize=12)
ax.set_ylabel("Proportion of Analyzable Genes", fontsize=12)
ax.set_title(
    "Gene sex bias\n"
    "(Purple = both-sex biased, Red = Female biased, Blue = Male biased)\n"
    "Left: per-species proportions | Right: conserved sex bias across species pairs/trio",
    fontsize=12,
)
ax.set_ylim(0, 1.0)
ax.grid(axis="y", alpha=0.3, linestyle="--")

# Legend
ax.legend(
    handles=[
        mpatches.Patch(facecolor=sexbias_color["B"], edgecolor="black", label="Both (B)"),
        mpatches.Patch(facecolor=sexbias_color["F"], edgecolor="black", label="Female (F)"),
        mpatches.Patch(facecolor=sexbias_color["M"], edgecolor="black", label="Male (M)"),
    ],
    loc="upper right",
    fontsize=10,
)

# Save PDF
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure4_panelC_02mdg.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.close()