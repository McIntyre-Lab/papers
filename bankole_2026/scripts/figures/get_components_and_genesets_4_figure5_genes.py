#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 09:40:02 2026

@author: mgaran
"""
import pandas as pd
PROJ= '/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/'
df = pd.read_csv(f"{PROJ}/zenodo/fiveSpecies_dmel6_full_annotation.csv")

gene_list = ["FBgn0030245","FBgn0000064","FBgn0017558","FBgn0015221","FBgn0015222","FBgn0026084"]

all_geneset = []
for gene in gene_list:
    geneset = df[df['geneID'] == gene]['genesetid'].iloc[0]
    all_geneset = all_geneset + [geneset]
    components = df[df['genesetid'] == geneset]['component_id'].astype(int).unique().tolist()
    comp_w_data = df[(df['genesetid'] == geneset) & (df['flag_jxnHash_in_data'] == 1)]['component_id'].astype(int).unique().tolist()
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print(gene)
    print(geneset)
    print(components)
    print(comp_w_data)
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("all genesets") 
print(all_geneset)

