#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 10:30:57 2026

@author: ammorse


add genesetID and componentID to datafile jxnHash

~/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files

component_map_by_node.csv
    jxnHash geneID source component_id
    ** will not have componentIDs for UJC not in annotation file
    
nodes_with_geneset.csv
    nodeid (geneid) genesetid
    ** should not be any missing    
    
"""


import pandas as pd

proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

comp = pd.read_csv(f"{proj}/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv")
comp.drop(columns="geneID", inplace=True)
#print(comp.columns)
    
geneset = pd.read_csv(f"{proj}/submission/supplementary_files/fiveSpecies_annotations/link_files/nodes_with_geneset.csv")
#print(geneset.columns)
geneset = geneset.rename(columns={'nodeid': 'geneID'})


species_list = ['dmel6', 'dsim2', 'dser1', 'dsan1', 'dyak2']
#species_list = ['dmel6']

for species in species_list:
    
    comp_species = comp.copy()
    comp_species = comp_species.rename(
        columns={'jxnhash': f"{species}_jxnHash"}
    ).reset_index(drop=True)

#    print(comp_species.columns)
    
    # import datafile
    datafile = pd.read_csv(f"{proj}/submission/supplementary_files/datafiles/datafile_jxnHash_{species}_almost.csv", low_memory=False)
    print(datafile.columns)
    
    ## merge datafile to geneset
    add_geneset = pd.merge(
        geneset,
        datafile, 
        on="geneID",
        how='outer',
        indicator='merge'
    )
    print(add_geneset['merge'].value_counts(dropna=False).sort_index())
#        left_only      127932
#        right_only          0
#        both          2286928 
    
    add_geneset = add_geneset[add_geneset['merge'].isin(['both'])]
    add_geneset.drop(columns="merge", inplace=True)
#    print(add_geneset.columns)

    # merge in component
    add_comp = pd.merge(
        comp_species,
        add_geneset, 
        on=f"{species}_jxnHash",
        how='outer',
        indicator='merge'
    )
    print(add_comp['merge'].value_counts(dropna=False).sort_index())
#       left_only      356930  only component
#       right_only    2264256  only in data
#       both            22672  in both

    #keep right and both, drop indicator col
    add_comp = add_comp[add_comp['merge'].isin(['right_only', 'both'])]
    add_comp.drop(columns="merge", inplace=True)

#    print(add_comp.columns)
    add_comp["geneID"].nunique()
    
    #** will not have componentIDs for UJC not in annotation file
    check = pd.crosstab(add_comp['component_id'], 
                      add_comp['flag_jxnHash_fiveSpecies_full_anno'], 
                      dropna=False, 
                      margins=True, 
                      margins_name="total").to_string()
          
    add_comp.to_csv(f"{proj}/zenodo/datafiles/datafile_jxnHash_{species}.csv", index=False)
    
    