#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
NETWORK = f"{PROJ}/zenodo/fiveSpecies_network_files"
OUT = f"{PROJ}/Figures_amm"

comp_df = pd.read_csv(f"{NETWORK}/component_info.csv")

species_colors = {
    'mel': '#966729',
    'sim': '#3F78C1',
    'san': '#28827A',
    'yak': '#717273',
    'ser': '#825CA6'
}

mask = (
    (comp_df['flag_expressed'] == 1) &
    (comp_df['comp_anno_in'] == '1_all_5_species')
)

expressed_in_cnts = comp_df.loc[mask, 'expressed_in'].value_counts()

order = [
    '4_mel_only',
    '4_sim_only',
    '4_san_only',
    '4_yak_only',
    '4_zser_only',
    '3_mel_sim',
    '3_yak_san',
    '6_other',
    '2_the4',
    '1_all_5'
]

expressed_in_cnts = expressed_in_cnts.reindex(order, fill_value=0)
expressed_in_cnts = expressed_in_cnts[expressed_in_cnts > 0]
        
label_map = {
    '4_mel_only': 'mel only',
    '4_sim_only': 'sim only',
    '4_san_only': 'san only',
    '4_yak_only': 'yak only',
    '4_zser_only': 'ser only',
    '3_mel_sim': 'mel + sim',
    '3_yak_san': 'yak + san',
    '6_other': 'other',
    '2_the4': 'mel,sim,san,yak',
    '1_all_5': 'all five'
}

# print if category is missing
missing = [cat for cat in order if cat not in expressed_in_cnts.index]
if missing:
    print("\nWARNING: The following categories were not found in the data and will not be plotted:")
    for cat in missing:
        print(f"  - {label_map.get(cat, cat)}")
        
        
style_map = {
    '4_mel_only':  (species_colors['mel'], 'white', ''),
    '4_sim_only':  (species_colors['sim'], 'white', ''),
    '4_san_only':  (species_colors['san'], 'white', ''),
    '4_yak_only':  (species_colors['yak'], 'white', ''),
    '4_zser_only': (species_colors['ser'], 'white', ''),
    '3_mel_sim':   (species_colors['mel'], species_colors['sim'], '//////'),
    '3_yak_san':   (species_colors['san'], species_colors['yak'], '//////'),
    '6_other':     ('#42A5F5', 'white', ''),
    '2_the4':      ('#1E88E5', 'white', ''),
    '1_all_5':     ('#0D47A1', 'white', '')
}

colors = [style_map[x][0] for x in expressed_in_cnts.index]
edgecolors = [style_map[x][1] for x in expressed_in_cnts.index]
hatches = [style_map[x][2] for x in expressed_in_cnts.index]

labels = [
    f"{label_map[idx]}\n({val / expressed_in_cnts.sum() * 100:.1f}%)"
    for idx, val in zip(expressed_in_cnts.index, expressed_in_cnts.values)
]

print("\nCOUNTS")
print(expressed_in_cnts)

print("\nPERCENT")
print((expressed_in_cnts / expressed_in_cnts.sum() * 100).round(3))

print("\nTOTAL")
print(expressed_in_cnts.sum())

fig, ax = plt.subplots(figsize=(10, 8))

wedges, _ = ax.pie(
    expressed_in_cnts.values,
    labels=labels,
    startangle=230,
    counterclock=False,
    colors=colors
)

for wedge, hatch, edgecolor in zip(wedges, hatches, edgecolors):
    wedge.set_hatch(hatch)
    wedge.set_edgecolor(edgecolor)
    wedge.set_linewidth(1.5 if hatch else 0.0) # linewidth of 0 makes hatches invisible

plt.title('Expressed components annotated in all 5 species by expression profile')
plt.axis('equal')

plt.savefig(
    f"{OUT}/figure2_panelB_02mdg.pdf",
    format='pdf',
    bbox_inches='tight',
    dpi=600
)
plt.close()

print(f"\nSaved: {OUT}/figure2_panelB_02mdg.pdf")