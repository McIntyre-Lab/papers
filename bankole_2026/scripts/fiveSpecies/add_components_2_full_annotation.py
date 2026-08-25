#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 11:02:27 2025

@author: ammorse

add components to full annotation files
use source col in component to link
"""

import pandas as pd
import gc 
gc.collect()

# Define file paths and species variable
proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

## note using dsim202 not dsimWXD
species_dct = {"dmel650":"dmel6", "dsan11":"dsan1", "dser11": "dser1", "dsim202": "dsim2", "dyak21": "dyak2"}
#species_dct = {"dmel650":"dmel6"}

for species_full,species in species_dct.items():
    
    # import component
    component = pd.read_csv(f"{proj}/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv", low_memory=False)
    #print(component.columns)
    #print(component["source"].value_counts())
    # rename jxnHash to {species}_jxnHash and merge to each anno
    component = component.rename(columns={"jxnhash": f"{species}_jxnHash"})
    
    # import anno file
    anno = pd.read_csv(f"{proj}/anno_testing/fiveSpecies_{species}_full_annotation_w_dataFlag_gffCompMM.csv", low_memory=False)
    #print(anno.columns)
    
    add_comp = pd.merge(
        component,
        anno,
        on = f"{species}_jxnHash",
        how = 'outer',
        indicator = 'merge')
    print(f"merge check for {species}")
    print(add_comp["merge"].value_counts())
        #left_only     302579
        #both           77023
        #right_only         0

    ## keep if in anno only and both 
    add_comp = add_comp[add_comp["merge"].isin(["right_only", "both"])].copy()
    print(f"row count for {species}")
    print(len(add_comp))    # 77023
    # drop merge col
    add_comp = add_comp.drop(columns=["merge"])
    
    # find any duplicate cols
    cols_xy = add_comp.columns[add_comp.columns.str.endswith(('_x', '_y'))]
    print(cols_xy)
    
    # check that geneID from both anno and component match
    print((add_comp['geneID_x'] == add_comp['geneID_y']).all())  
    ## True so drop 1 and rename other
    add_comp = add_comp.rename(columns={"geneID_x": "geneID"})
    add_comp = add_comp.drop(columns=["geneID_y"])

    # check that cat_dsan1_transcriptID from both anno and component match
    add_comp[f"cat_{species}_transcriptID_x"].equals(add_comp[f"cat_{species}_transcriptID_y"])
    ## True so rename and drop on1    
    add_comp = add_comp.rename(columns={f"cat_{species}_transcriptID_x": f"cat_{species}_transcriptID"})
    add_comp = add_comp.drop(columns=[f"cat_{species}_transcriptID_y"])
       
    # save to csv
    add_comp.to_csv(f"{proj}/anno_testing/fiveSpecies_{species}_full_annotation_w_component.csv", index=False)    
    