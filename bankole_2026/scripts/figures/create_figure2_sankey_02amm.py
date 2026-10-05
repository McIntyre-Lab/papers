#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import plotly.graph_objects as go
from collections import defaultdict

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
TABLES = f"{PROJ}/Tables"
OUT = f"{PROJ}/Figures_amm"

sankey_df = pd.read_csv(f"{TABLES}/cnt_origin_anno_xpress_byspecies.csv")
sankey_df['flag_panelA_sankey'] = (sankey_df['expressed_in'] != '7_not_expressed').astype(int)
sankey_df["expressed_in"].isna().any()

# categories to keep, collapse others into 'other'
keep = ['1_all_5_species', '2_the4', '4_zseronly']
sankey_df['comp_anno_in'] = sankey_df['comp_anno_in'].where(sankey_df['comp_anno_in'].isin(keep), 'other')

comp_anno_flags = sankey_df.groupby('comp_anno_in')['flag_panelA_sankey'].max().to_dict()
expressed_flags = sankey_df.groupby('expressed_in')['flag_panelA_sankey'].max().to_dict()
## keep where at least 1 rows is expressed 
comp_annos = [cat for cat, flag in comp_anno_flags.items() if flag == 1]
expressed_vals = [cat for cat, flag in expressed_flags.items() if flag == 1]

comp_anno_labels = [f"Comp Anno: {x}" for x in comp_annos]
expressed_labels = [f"Expressed: {x}" for x in expressed_vals]
all_labels = comp_anno_labels + expressed_labels
label_to_idx = {label: idx for idx, label in enumerate(all_labels)}

link_counts = defaultdict(int)
for _, row in sankey_df.iterrows():
    if row['flag_panelA_sankey'] == 1:
        source_idx = label_to_idx[f"Comp Anno: {row['comp_anno_in']}"]
        target_idx = label_to_idx[f"Expressed: {row['expressed_in']}"]
        link_counts[(source_idx, target_idx)] += row['COUNT']

sources = [k[0] for k in link_counts]
targets = [k[1] for k in link_counts]
values = list(link_counts.values())

fig = go.Figure(data=[go.Sankey(
    node=dict(pad=15, thickness=20, line=dict(color='black', width=0.5), label=all_labels, color='blue'),
    link=dict(source=sources, target=targets, value=values)
)])
fig.update_layout(title="Sankey Plot: Comp Anno → Expressed", font=dict(size=10), height=800, width=1200)
outfile = f"{OUT}/figure2_sankey_02mdg.pdf"
fig.write_image(outfile, format='pdf')
print(f"Saved: {outfile}")
