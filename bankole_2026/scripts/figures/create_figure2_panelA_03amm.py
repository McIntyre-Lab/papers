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

category_map = {
    '4_yakonly': 'one_species_yak/mel/san/sim',
    '4_melonly': 'one_species_yak/mel/san/sim',
    '4_sanonly': 'one_species_yak/mel/san/sim',
    '4_simonly': 'one_species_yak/mel/san/sim',
    '4_zseronly': 'one_species_ser',
    '3_mel_sim': 'other',
    '3_yak_san': 'other',
    '6_other': 'other',
    '2_the4': 'mel,sim,san,yak',
    '1_all_5_species': 'mel,sim,san,yak,ser'
}

expressed = comp_df['flag_expressed'] == 1

raw_cnts = comp_df.loc[expressed, 'comp_anno_in'].value_counts()

mapped_index = raw_cnts.index.map(lambda x: category_map.get(x, x))

combined_cnts = (
    pd.Series(index=mapped_index, data=raw_cnts.values)
    .groupby(level=0)
    .sum()
)

order = [
    'one_species_yak/mel/san/sim',
    'one_species_ser',
    'other',
    'mel,sim,san,yak',
    'mel,sim,san,yak,ser'
]

combined_cnts = combined_cnts.reindex(order, fill_value=0)

print("\nCOUNTS")
print(combined_cnts)

print("\nPERCENT")
print((combined_cnts / combined_cnts.sum() * 100).round(3))

print("\nTOTAL")
print(combined_cnts.sum())

plt.figure(figsize=(10, 8))
plt.pie(
    combined_cnts.values,
    labels=combined_cnts.index,
    autopct='%1.1f%%',
    startangle=210,
    counterclock=False,
    colors=['#E3F2FD', '#90CAF9', '#42A5F5', '#1E88E5', '#0D47A1']
)

plt.title('Expressed components annotated in (comp_anno_in)')
plt.axis('equal')

plt.savefig(
    f"{OUT}/figure2_panelA_02mdg.pdf",
    format='pdf',
    bbox_inches='tight',
    dpi=600
)
plt.close()

print(f"\nSaved: {OUT}/figure2_panelA_02mdg.pdf")