#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 16:39:22 2026

@author: ammorse


add genesetID to full anno 

"""


import pandas as pd

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

# note using dsim202 not dsimWXD
species_dct = {"dmel650":"dmel6", "dsan11":"dsan1", "dser11": "dser1", "dsim202": "dsim2", "dyak21": "dyak2"}
#species_dct = {"dmel650":"dmel6"}

for species_full,species in species_dct.items():
    
    # imp genesetid
    node = pd.read_csv(f"{PROJ}/submission/supplementary_files/fiveSpecies_annotations/link_files/nodes_with_geneset.csv", low_memory=False)
    print(node.columns)
    # rename nodeid to geneID
    node = node.rename(columns={"nodeid": "geneID"})
    
    # import anno file
    anno = pd.read_csv(f"{PROJ}/anno_testing/fiveSpecies_{species}_full_annotation_w_component.csv", low_memory=False)
    print(anno.columns)
    print(len(anno))
    
    add_geneset = pd.merge(
        node,
        anno,
        on = 'geneID',
        how = 'outer',
        indicator = 'merge')
    print(f"merge check for {species}")
    print(add_geneset["merge"].value_counts())
    #    merge check for dmel6
    #    merge
    #    left_only     125700  node only 
    #    both           75986  both
    #    right_only            anno only
    
    ## keep if in anno only and both 
    add_geneset = add_geneset[add_geneset["merge"].isin(["right_only", "both"])].copy()
    print(f"row count for {species}")
    print(len(add_geneset))    # 77023
    # drop merge col
    add_geneset = add_geneset.drop(columns=["merge"])
    
        # find any duplicate cols
    cols_xy = add_geneset.columns[add_geneset.columns.str.endswith(('_x', '_y'))]
    print(cols_xy)

      
    # save to csv
    add_geneset.to_csv(f"{PROJ}/anno_testing/fiveSpecies_{species}_full_annotation.csv", index=False)    
    
    
    