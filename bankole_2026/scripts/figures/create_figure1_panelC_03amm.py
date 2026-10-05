#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 18:47:55 2026

@author: mgaran
"""
import pandas as pd
import plotly.graph_objects as go
from collections import defaultdict

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
#PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
OUT = f"{PROJ}/Figures_amm"

# Read count file
df = pd.read_csv(f"{PROJ}/Tables/cnt_origin_anno_xpress_byspecies.csv")

# Print total COUNT grouped by origin
print("\nORIGIN COUNTS:")
origin_counts = df.groupby("origin")["COUNT"].sum().sort_values(ascending=False)
print(origin_counts.to_string())
print("TOTAL:", origin_counts.sum())

# Print total COUNT grouped by comp_anno_in
print("\nCOMP_ANNO_IN COUNTS:")
comp_anno_counts = df.groupby("comp_anno_in")["COUNT"].sum().sort_values(ascending=False)
print(comp_anno_counts.to_string())
print("TOTAL:", comp_anno_counts.sum())

# Create unique node labels
origins = df["origin"].unique().tolist()
comp_annos = df["comp_anno_in"].unique().tolist()

all_labels = (
    [f"Origin: {x}" for x in origins]
    + [f"Comp Anno: {x}" for x in comp_annos]
)

# Map node labels to index
label_to_idx = {
    label: idx
    for idx, label in enumerate(all_labels)
}

# Aggregate COUNT per origin + comp_anno_in pair
link_counts = defaultdict(int)

for _, row in df.iterrows():
    source = label_to_idx[f"Origin: {row['origin']}"]
    target = label_to_idx[f"Comp Anno: {row['comp_anno_in']}"]
    link_counts[(source, target)] += row["COUNT"]

sources = [x[0] for x in link_counts]
targets = [x[1] for x in link_counts]
values = list(link_counts.values())

# Plot Sankey
fig = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color="blue",
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
            ),
        )
    ]
)

fig.update_layout(
    title="Sankey Plot: Origin -> Comp Anno",
    font=dict(size=12),
    height=800,
    width=1200,
)

fig.write_image(
    f"{OUT}/figure1_panelC_02mdg.pdf",
    format="pdf",
    scale=4,
)