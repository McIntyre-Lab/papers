#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 17:34:11 2025

@author: nkeil
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

# === I/O dirs ===
ind  = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis"
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis"

# File produced by script #1 (contains counts, total, and perc_* columns)
table_csv = os.path.join(ind, "structural_category_plus_counts_by_sample.csv")

# Stacked bar order (bottom -> top)
#plot_order = ["FSM", "NIC", "NNC_annoERP", "ISM", "NNC_no_annoERP"]
plot_order = ["FSM", "NIC", "NNC_annoERP", "ISM"]


# Colors (NNC color reused for both sub-categories)
category_color_palette = {
    "FSM": "#6BAED6",
    "ISM": "#FC8D59",
    "NIC": "#78C679",
    "NNC": "#EE6A50",
    "GENIC": "#969696",
    "AS": "#66C2A4",
    "FUS": "#ffc125",
    "INTER": "#e9967a",
    "GI": "#41B6C4",
}
color_map = {
    "FSM": category_color_palette["FSM"],
    "NIC": category_color_palette["NIC"],
    "NNC_annoERP": category_color_palette["NNC"],
    "ISM": category_color_palette["ISM"],
    "NNC_no_annoERP": category_color_palette["NNC"],
}

# Hatching: bottom 3 = both diagonals ('x'), ISM = '/', NNC_no_annoERP = none
hatch_map = {
    "FSM": "",
    "NIC": "",
    "NNC_annoERP": "",
    "ISM": "",
    "NNC_no_annoERP": "o",
}

# --- load table ---
df = pd.read_csv(table_csv, index_col=0)

# Ensure needed percentage columns exist (fallback to zero if missing)
for cat in plot_order:
    col = f"perc_{cat}"
    if col not in df.columns:
        df[col] = 0.0

df_perc = df[[f"perc_{cat}" for cat in plot_order]].copy()

# --- plot ---
samples = df_perc.index.tolist()
x = range(len(samples))

fig_w = max(8, len(samples) * 0.5)
fig, ax = plt.subplots(figsize=(fig_w, 5))

bottom = [0.0] * len(samples)
for cat in plot_order:
    heights = df_perc[f"perc_{cat}"].values
    bars = ax.bar(
        x, heights, bottom=bottom,
        label=cat,
        color=color_map.get(cat, "#cccccc"),
        edgecolor="black",
        linewidth=0.4,
    )
    hatch = hatch_map.get(cat, "")
    if hatch:
        for b in bars:
            b.set_hatch(hatch)
    bottom = [b + h for b, h in zip(bottom, heights)]

ax.set_ylabel("Percent of reads")
ax.set_ylim(0, 100)
ax.set_xticks(list(x))
ax.set_xticklabels(samples, rotation=75, ha="right")
ax.legend(title="Category", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
ax.set_title("Per-sample structural category plus (%)")

plt.tight_layout()
plot_out = os.path.join(outd, "struct_category_plus_percent_stacked_02ndk.png")
plt.savefig(plot_out, dpi=300)
plt.close()
#plt.show()

print(f"✓ Wrote plot: {plot_out}")
