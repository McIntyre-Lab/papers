#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
exclude ISM ujc
exclude ujc < 10 reads
only multiexon genes


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

category_order = ["1", "2", "3", "4", "5", "6+"]
tier_order = ["anno_UJC", "anno_ERP", "novel_ERP"]
color_by_tier = {"anno_UJC": dark, "anno_ERP": medium, "novel_ERP": light}

panel_b_data = {}
num_gene_ge2 = {}

for species in species_list:
    df = pd.read_csv(
        f"{SUMMARY}/gene_summary_min10reads_{species}.csv",
        usecols=[
            "geneID", "numExon_GM", "num_ujc", "num_ism_ujc",
            "num_anno_ujc", "num_ujc_w_anno_erp", "num_ujc_w_novel_erp"
        ],
        low_memory=False
    )

    df = df[df["numExon_GM"] >= 2].copy()
    num_gene_ge2[species] = len(df)

    df["num_nonISM_UJC"] = df["num_ujc"] - df["num_ism_ujc"]
    df = df[df["num_nonISM_UJC"] >= 1].copy()

    df["num_nonISM_bin"] = df["num_nonISM_UJC"].apply(
        lambda x: "6+" if x >= 6 else str(int(x))
    )
    
    
    ## tier assignment with sequential overwriting - last assignment wins
        # anno_ujc is last and will overwrite anno_ERP and novel_ERP....    
    df["annotation_tier"] = ""
    df.loc[df["num_ujc_w_novel_erp"] >= 1, "annotation_tier"] = "novel_ERP"
    df.loc[df["num_ujc_w_anno_erp"] >= 1, "annotation_tier"] = "anno_ERP"
    df.loc[df["num_anno_ujc"] >= 1, "annotation_tier"] = "anno_UJC"

    # check for genes with all 0
    unassigned = (df["annotation_tier"] == "").sum()
    if unassigned > 0:
        print(f"  WARNING: {unassigned} genes in {species} have no annotation tier assigned")
        
    panel_b_data[species] = {}
    for category in category_order:
        df_cat = df[df["num_nonISM_bin"] == category]
        panel_b_data[species][category] = {
            tier: int((df_cat["annotation_tier"] == tier).sum())
            for tier in tier_order
        }

print("\nGENES GOING INTO PLOT PER SPECIES")
for species in species_list:
    total = sum(
        panel_b_data[species][cat][tier]
        for cat in category_order
        for tier in tier_order
    )
    print(f"  {species}: {total:,}")
    
fig, ax = plt.subplots(figsize=(16, 8))
x = np.arange(len(category_order))
bar_width = 0.12
offset_positions = np.arange(len(species_list)) - (len(species_list) - 1) / 2

for i, species in enumerate(species_list):
    x_pos = x + offset_positions[i] * bar_width
    bottom = np.zeros(len(category_order))

    for tier in tier_order:
        counts = [panel_b_data[species][category][tier] for category in category_order]
        ax.bar(
            x_pos, counts, bar_width, bottom=bottom,
            color=color_by_tier[tier][species],
            edgecolor="black", linewidth=0.5
        )
        bottom += np.array(counts)

ax.set_xlabel("Number of non-ISM UJCs per gene", fontsize=12)
ax.set_ylabel("Number of genes", fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(category_order, fontsize=11)
ax.tick_params(axis="y", labelsize=11)
ax.grid(axis="y", alpha=0.3, linestyle="--")

plt.tight_layout()
plt.savefig(f"{OUT}/figure3_panelB_02mdg.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.close()

print("numGene_w_numExon_GM_ge2\tdmel\tdsim\tdyak\tdsan\tser")
print("count\t" + "\t".join(str(num_gene_ge2[s]) for s in species_list))
