#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  6 20:28:48 2025

@author: nkeil
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

# Set paths to directories for sqanti reads vs 5 species
rmg_5sp_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_facet_sex"
rlr_5sp_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_facet_sex"
axk_5sp_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_facet_sex"

# Set paths to directories for sqanti reads vs "native" annotations
rmg_nat_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rmg_dros_data/sqantiReads_vs_native_rmg_facet_sex"
rlr_nat_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_rlr_head_data/sqantiReads_vs_native_facet_sex"
axk_nat_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/lmm_axk_head_data/sqantiReads_vs_native_kopp_facet_sex"


#Set paths to structural category files
mel_nat_file = os.path.join(rmg_nat_dir, "rmg_mel_comb_TRs_facet_sex_structural_category_percs.csv")
sim_nat_file = os.path.join(rmg_nat_dir, "rmg_sim_comb_TRs_facet_sex_structural_category_percs.csv")
yak_nat_file = os.path.join(rlr_nat_dir, "TO_dyak_sex_structural_category_percs.csv")
san_nat_file = os.path.join(rlr_nat_dir, "TO_dsan_sex_structural_category_percs.csv")
ser_nat_file = os.path.join(axk_nat_dir, "ser_facet_sex1_structural_category_percs.csv")

mel_5sp_file = os.path.join(rmg_5sp_dir, "rmg_mel_comb_TRs_facet_sex_structural_category_percs.csv")
sim_5sp_file = os.path.join(rmg_5sp_dir, "rmg_sim_comb_TRs_facet_sex_structural_category_percs.csv")
yak_5sp_file = os.path.join(rlr_5sp_dir, "TO_dyak_comb_TRs_facet_sex_structural_category_percs.csv")
san_5sp_file = os.path.join(rlr_5sp_dir, "TO_dsan_comb_TRs_facet_sex_structural_category_percs.csv")
ser_5sp_file = os.path.join(axk_5sp_dir, "kopp_dser_comb_TRs_facet_sex_structural_category_percs.csv")

#Set up file directory
species_files = {
    "mel": {
        "nat_file": mel_nat_file,
        "5sp_file": mel_5sp_file,
        "tag": "dmel6"
    },
    "sim": {
        "nat_file": sim_nat_file,
        "5sp_file": sim_5sp_file,
        "tag": "dsim2"
    },
    "yak": {
        "nat_file": yak_nat_file,
        "5sp_file": yak_5sp_file,
        "tag": "dyak2"
    },
    "san": {
        "nat_file": san_nat_file,
        "5sp_file": san_5sp_file,
        "tag": "dsan1"
    },
    "ser": {
        "nat_file": ser_nat_file,
        "5sp_file": ser_5sp_file,
        "tag": "dser1"
    }
}

#Set up list
combined_data = []

#Loop through 5 species

for species, files in species_files.items():
    file_nat = files["nat_file"]
    file_5sp = files["5sp_file"]

    df_nat = pd.read_csv(file_nat)
    df_5sp = pd.read_csv(file_5sp)

    # Check column alignment
    assert list(df_nat.columns) == list(df_5sp.columns), f"Column mismatch in {species}"

    sort_cols = ['sampleID']
    df_nat = df_nat.sort_values(by=sort_cols).reset_index(drop=True)
    df_5sp = df_5sp.sort_values(by=sort_cols).reset_index(drop=True)

    assert df_nat[sort_cols].equals(df_5sp[sort_cols]), f"Row mismatch in {species}"

    # Create combined dataframe with suffixes
    numeric_cols = df_nat.select_dtypes(include='number').columns
    df_comb = pd.DataFrame()
    df_comb[sort_cols] = df_nat[sort_cols]
    
    for col in numeric_cols:
        df_comb[f"{col}_native"] = df_nat[col]
        df_comb[f"{col}_5species"] = df_5sp[col]
        df_comb[f"{col}_diff"] = df_5sp[col] - df_nat[col]

    df_comb["species"] = species
    df_comb["tag"] = files["tag"]
    
    combined_data.append(df_comb)

# Combine all species into one DataFrame
df_combined = pd.concat(combined_data, ignore_index=True)

## Create new combined GENIC_INTER columns
df_combined["GENIC_INTER_native"] = df_combined["GENIC_native"] + df_combined["INTER_native"]
df_combined["GENIC_INTER_5species"] = df_combined["GENIC_5species"] + df_combined["INTER_5species"]
df_combined["GENIC_INTER_diff"] = df_combined["GENIC_INTER_5species"] - df_combined["GENIC_INTER_native"]

## Create new combined FSM+ISM+NIC+NNC columns
df_combined["FSM_ISM_NIC_NNC_native"] = df_combined["FSM_native"] + df_combined["ISM_native"] + df_combined["NIC_native"] + df_combined["NNC_native"]
df_combined["FSM_ISM_NIC_NNC_5species"] =  df_combined["FSM_5species"] + df_combined["ISM_5species"] + df_combined["NIC_5species"] + df_combined["NNC_5species"]
df_combined["FSM_ISM_NIC_NNC_diff"] = df_combined["FSM_ISM_NIC_NNC_5species"] - df_combined["FSM_ISM_NIC_NNC_native"]

# Reorder columns
first_cols = ['sampleID', 'species', 'tag']
other_cols = [col for col in df_combined.columns if col not in first_cols]
df_combined = df_combined[first_cols + other_cols]

# Save combined dataframe with native, 5species, and diff
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables"

#Categories of intterest
plot_categories = ['FSM', 'ISM','NIC','NNC' ,'GENIC', 'INTER','GENIC_INTER','FSM_ISM_NIC_NNC']

# Limit to relevant columns only
base_cols = ['sampleID', 'species', 'tag']
suffixes = ['_native', '_5species', '_diff']
selected_cols = [f"{cat}{suf}" for cat in plot_categories for suf in suffixes]
df_export = df_combined[base_cols + selected_cols]
df_export.to_csv(os.path.join(outd,"sqanti_compare_native_vs_5species.csv"))

### PLOTTING SECTION ###

#Barplot

plot_diff_cols = [f"{cat}_diff" for cat in plot_categories]

# Set up plot grid
n_cols = 2
n_plots = len(plot_diff_cols)
n_rows = math.ceil(n_plots / n_cols)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 4 * n_rows))
axes = axes.flatten()

# Plot each diff category
for idx, col in enumerate(plot_diff_cols):
    ax = axes[idx]
    
    # Bar plot
    ax.bar(df_combined['sampleID'], df_combined[col])
    ax.set_title(col.replace("_diff", ""))
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Difference")
    ax.set_xticks(range(len(df_combined['sampleID'])))
    ax.set_xticklabels(df_combined['sampleID'], rotation=45, ha='right')

# Turn off unused axes
for idx in range(n_plots, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.savefig("/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_native_vs_5species_bar.png")
plt.show()

#SCATTERPLOT

# Prepare long-form DataFrame for scatterplot
scatter_data = []

for cat in plot_categories:
    for _, row in df_combined.iterrows():
        sex = "F" if "F" in row["sampleID"] else "M"
        scatter_data.append({
            "category": cat,
            "native": row[f"{cat}_native"],
            "5species": row[f"{cat}_5species"],
            "species": row["species"],
            "sex": sex,
            "sampleID": row["sampleID"]
        })

df_scatter = pd.DataFrame(scatter_data)

# Define colors for species
species_palette = {
    "mel": "#44a043",
    "sim": "#3f78c1",
    "yak": "#7c7c7c",
    "san": "#28827a",
    "ser": "#8252a6"
}

rep_markers = {
    "1-3": "o",  # reps 1–3
    "4-6": "s"    # reps 4–6
}

# Set up plot grid
fig, axes = plt.subplots(n_cols, n_rows, figsize=(14, 10))
axes = axes.flatten()

# Track seen labels to avoid duplicates in legend
seen_labels = set()

for idx, category in enumerate(plot_categories):
    ax = axes[idx]

    for _, row in df_combined.iterrows():
        sample_id = row["sampleID"]
        species = row["species"]

        # Extract replicate number (1–6), default to 1 if not found
        rep_num = next((int(c) for c in sample_id if c.isdigit()), 1)

        # Determine rep group for labeling/marker shape
        if species in ["yak", "san"]:
            rep_label = f"{species} {rep_num}"
            marker_key = "1-3"  # still use "o" as marker
        else:
            rep_group = "1-3" if rep_num <= 3 else "4-6"
            rep_label = f"{species} {rep_group}"
            marker_key = rep_group

        # Only add label once to avoid duplicates
        show_label = rep_label not in seen_labels
        if show_label:
            seen_labels.add(rep_label)

        ax.scatter(
            row[f"{category}_native"],
            row[f"{category}_5species"],
            color=species_palette.get(species, "gray"),
            marker=rep_markers.get(marker_key, "x"),
            label=rep_label if show_label else None,
            s=60,
            edgecolor="black"
        )

    ax.set_title(f"{category}:  5-Species vs Native")
    ax.set_xlabel("Native %")
    ax.set_ylabel("5-Species %")

    # 45-degree reference line
    max_val = max(df_combined[f"{category}_native"].max(), df_combined[f"{category}_5species"].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3)

# Add one global legend (from the first axis)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=4, fontsize="medium")

# Hide unused panels
if len(plot_categories) < len(axes):
    for j in range(len(plot_categories), len(axes)):
        axes[j].axis('off')
        
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_native_vs_5species_scatte.png", bbox_inches='tight')
plt.show()
        
#New plot FSM_ISM_NIC_NNC_diff vs GENIC_INTERGENIC_diff
fig, ax = plt.subplots(figsize=(8, 6))

# reset seen labels for this figure
seen_labels = set()

for _, row in df_combined.iterrows():
    sample_id = row["sampleID"]
    species = row["species"]

    rep_num = next((int(c) for c in sample_id if c.isdigit()), 1)

    if species in ["yak", "san"]:
        rep_label = f"{species} {rep_num}"
        marker_key = "1-3"
    else:
        rep_group = "1-3" if rep_num <= 3 else "4-6"
        rep_label = f"{species} {rep_group}"
        marker_key = rep_group

    show_label = rep_label not in seen_labels
    if show_label:
        seen_labels.add(rep_label)

    ax.scatter(
        row["FSM_ISM_NIC_NNC_diff"],
        row["GENIC_INTER_diff"],
        color=species_palette.get(species, "gray"),
        marker=rep_markers.get(marker_key, "x"),
        label=rep_label if show_label else None,   # only add new labels
        s=70,
        edgecolor="black"
    )

ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

ax.set_title("FSM+ISM+NIC+NNC (Diff) vs GENIC+INTER (Diff)")
ax.set_xlabel("FSM+ISM+NIC+NNC Difference (5sp - native)")
ax.set_ylabel("GENIC+INTER Difference (5sp - native)")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc="upper left", fontsize="small")

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis/FSM_ISM_NIC_NNC_diff_vs_GENIC_INTERGENIC_diff.png", bbox_inches='tight')
plt.show()

# Barplot - FSM, NIC, GENIC-INTER on same plot
import numpy as np

bar_categories = ['FSM', 'NIC', 'GENIC_INTER']
plot_cols = [f"{cat}_diff" for cat in bar_categories]
colors = ['#6BAED6', '#78C679', '#969696']  # SQANTI3 colors
labels = ['FSM', 'NIC', 'GENIC_INTER']

x = np.arange(len(df_combined['sampleID']))
bar_width = 0.25
offsets = np.linspace(-bar_width, bar_width, len(plot_cols))

fig, ax = plt.subplots(figsize=(max(10, len(x) * 0.6), 6))

# Create grouped bars
for i, (col, color, label, offset) in enumerate(zip(plot_cols, colors, labels, offsets)):
    ax.bar(x + offset, df_combined[col], width=bar_width,
           label=label, color=color, edgecolor='black', linewidth=0.7)

# Styling
ax.set_title("FSM, NIC, and GENIC_INTER Differences (5 species - native)")
ax.set_xlabel("Sample ID")
ax.set_ylabel("Difference ( 5 species - Native (%) )")
ax.set_xticks(x)
ax.set_xticklabels(df_combined['sampleID'], rotation=90, ha='center')
ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax.legend(frameon=False, ncol=3, loc='upper right')

plt.tight_layout()
plt.savefig("/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_native_vs_5species_diff_bar_sameplot.png",
            dpi=300, bbox_inches='tight')

# keep text editable (not outlines)
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["text.usetex"] = False

plt.savefig(
    "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_native_vs_5species_diff_bar_sameplot.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
)

plt.show()

