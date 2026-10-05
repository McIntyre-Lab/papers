#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""


import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BASE_DIR = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures_amm"
SUMMARY_DIR = f"{BASE_DIR}/zenodo/summary_files"

sex_bias_species_list = ["dmel", "dsim", "dser"]

dark = {"dmel": "#966729", "dsim": "#3F78C1", "dsan": "#28827A", "dyak": "#717273", "dser": "#825CA6"}
medium = {"dmel": "#D4B08C", "dsim": "#9FC1E8", "dsan": "#8FCCC7", "dyak": "#B8B9BA", "dser": "#C3A9D3"}
light = {"dmel": "#F4D8A8", "dsim": "#C8E2F8", "dsan": "#C6E9E5", "dyak": "#D9DADB", "dser": "#E1D8F0"}

gene_model_size_order = ["short", "med", "long"]

sexbias_stacks = {
    "B": {"label": "B", "palette": dark},
    "F_or_M": {"label": "F or M", "palette": medium},
    "unbiased": {"label": "Unbiased", "palette": light},
}

bar_positions = {
    "dmel": {"short": 0, "med": 1, "long": 2},
    "dsim": {"short": 4, "med": 5, "long": 6},
    "dser": {"short": 8, "med": 9, "long": 10},
}

panel_a_data = {}

for species in sex_bias_species_list:
    df = pd.read_csv(
        f"{SUMMARY_DIR}/gene_summary_{species}.csv",
        usecols=["geneID", "sexBias", "numExon_GM", "num_ujc_analyzable"],
        low_memory=False,
    )
    
    ## make sure numExon_GM == 0 (multiexon) or NaN not included
    ## including monoexons!
    df = df[df["numExon_GM"].notna() & (df["numExon_GM"] > 0)]
    
    df["gm_size"] = ""
    df.loc[df["numExon_GM"] <= 3, "gm_size"] = "short"
    df.loc[(df["numExon_GM"] > 3) & (df["numExon_GM"] <= 9), "gm_size"] = "med"
    df.loc[df["numExon_GM"] > 9, "gm_size"] = "long"

    print(df["sexBias"].unique())
    
    df["sexBias_group"] = ""
    df.loc[df["sexBias"] == "B", "sexBias_group"] = "B"
    df.loc[df["sexBias"].isin(["F", "M"]), "sexBias_group"] = "F_or_M"
    df.loc[df["sexBias"] == "unbiased", "sexBias_group"] = "unbiased"

    analyzable = df["num_ujc_analyzable"] >= 1
    
                        
    panel_a_data[species] = {}

    for gm_size in gene_model_size_order:
        denominator = analyzable & (df["gm_size"] == gm_size)
        denominator_count = denominator.sum()
        panel_a_data[species][gm_size] = {}

        for sexbias_group in sexbias_stacks:
            count = (denominator & (df["sexBias_group"] == sexbias_group)).sum()
            panel_a_data[species][gm_size][sexbias_group] = (
                count / denominator_count if denominator_count > 0 else 0
            )

    print(f"\n=== {species} counts ===")
    print(f"  Raw rows in CSV:                        {len(pd.read_csv(f'{SUMMARY_DIR}/gene_summary_{species}.csv', usecols=['geneID','sexBias','numExon_GM','num_ujc_analyzable'], low_memory=False))}")
    print(f"  After numExon_GM > 0 and notna:         {len(df)}")
    
    print(f"  With valid sexBias_group:                {(df['sexBias_group'] != '').sum()}")
    print(f"  Analyzable (num_ujc_analyzable >= 1):    {(df['num_ujc_analyzable'] >= 1).sum()}")
    print(f"  Analyzable AND valid sexBias_group:      {((df['num_ujc_analyzable'] >= 1) & (df['sexBias_group'] != '')).sum()}")


fig, ax = plt.subplots(figsize=(10, 6))

for species, positions in bar_positions.items():
    for gm_size, xpos in positions.items():
        bottom = 0
        for sexbias_group, settings in sexbias_stacks.items():
            proportion = panel_a_data[species][gm_size][sexbias_group]
            ax.bar(
                xpos, proportion, bottom=bottom,
                color=settings["palette"][species],
                width=0.8, edgecolor="black", linewidth=0.5,
            )
            bottom += proportion

ax.set_xticks([0, 1, 2, 4, 5, 6, 8, 9, 10])
ax.set_xticklabels([
    "short", "med", "long",
    "short", "med", "long",
    "short", "med", "long",
])
ax.set_ylim(0, 1)
ax.set_ylabel("Proportion of genes")
ax.set_xlabel("Gene model size")

ax.text(1, -0.14, "dmel", transform=ax.get_xaxis_transform(), ha="center", va="top")
ax.text(5, -0.14, "dsim", transform=ax.get_xaxis_transform(), ha="center", va="top")
ax.text(9, -0.14, "dser", transform=ax.get_xaxis_transform(), ha="center", va="top")

ax.axvline(3, color="black", linewidth=0.8)
ax.axvline(7, color="black", linewidth=0.8)

ax.legend(
    handles=[
        Patch(facecolor="#E1E1E1", edgecolor="black", label="Light: unbiased"),
        Patch(facecolor="#999999", edgecolor="black", label="Medium: F or M"),
        Patch(facecolor="#333333", edgecolor="black", label="Dark: B"),
    ],
    title="Sex bias",
)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure4_panelA_02mdg.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.close()
