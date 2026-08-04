#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 16:02:10 2025

@author: ammorse

add flag_gene_misMatch

mel  30508
san  25717
yak  27479
ser  21009
sim  21441

"""

import pandas as pd
import gc 
gc.collect()

# Define file paths and species variable
proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"  # Base file path

## note using dsim202 not dsimWXD
species_dct = {"dmel650":"dmel6", "dsan11":"dsan1", "dser11": "dser1", "dsim202": "dsim2", "dyak21": "dyak2"}
#species_dct = {"dsan11":"dsan1"}

for species_full,species in species_dct.items():
    
    # import anno file
    anno = pd.read_csv(f"{proj}/submission/supplementary/fiveSpecies_{species}_full_annotation_w_dataFlag.csv", low_memory=False)
    #print(anno.columns)
    print(f"rows in {species} anno: {len(anno)}")   # 75986 rows
    print(f"# uniq geneIDs in {species} anno:", anno["geneID"].nunique())
    num_dupes = anno[f"{species}_jxnHash"].duplicated().sum()
    print(f"# of dup jxnHash 4 {species}:", num_dupes)

        # # uniq geneIDs in dmel6 anno: 18929
    
    ## explode cat_dmel6_transcriptID
    anno_exp = anno.assign(
    transcriptID=anno[f'cat_{species}_transcriptID'].str.split('|')
    ).explode('transcriptID')
    print(f"rows in {species} anno exploded: {len(anno_exp)}")   # 79787 rows
    num_dupes = anno_exp["transcriptID"].duplicated().sum()
    print(f"# of dup transcriptID in anno_exp 4 {species}:", num_dupes)
        # 51305 for san
    
    # import gene flag file
    flag = pd.read_csv(f"{proj}/cross_species_link_files/{species_full}_2_{species}_compare_transcript_gene_assignment.csv", low_memory=False)
    #print(flag.columns)  
    print(f"Total rows in {species} flag file: {len(flag)}")
        # Total rows in dmel6 flag file: 34309
    print(f"Unique transcriptIDs: {flag['transcriptID'].nunique()}")
        # Unique transcriptIDs: 34309
    num_dupes = flag["transcriptID"].duplicated().sum()
    print(f"# of dup transcriptID in flag 4 {species}:", num_dupes)    
    
    # keep only cols need
    keep = ['transcriptID', 'flag_geneMismatch']  
    flag = flag[keep]
    flag = flag.drop_duplicates()
    num_dupes = flag["transcriptID"].duplicated().sum()
    print(f"# of dup transcriptID in flag 4 {species}:", num_dupes)
        # 0 dups in flag for dsan 
        
    # merge flag into anno_exp -- 1 to many merge....     
    anno_flag = pd.merge(
        flag,
        anno_exp,
        left_on = 'transcriptID',
        right_on = 'transcriptID',
        how = 'outer',
        indicator = 'merge')
    counts = anno_flag["merge"].value_counts()
    # Extract counts safely (0 if not present)
    right_only = counts.get("right_only", 0)
    left_only = counts.get("left_only", 0)
    both = counts.get("both", 0)
    print(f"L_only: {left_only}, R_only: {right_only}, B: {both}, R+B: {right_only + both}")
        # L_only: 0, R_only: 45478, B: 34309, R+B: 79787
        
    ##  keep if right only and both
    added_flag = anno_flag[anno_flag["merge"].isin(["right_only", "both"])].copy()

    ## collapse back to cat_{species}_transcriptID by {species}_jxnHash
    collapse = (
        added_flag
        .groupby(f"{species}_jxnHash", as_index=False)
        .agg({
            "transcriptID": lambda x: "|".join(x.dropna().astype(str)),
            "flag_geneMismatch": "max"   # take max flag per jxnHash
        })
        .rename(columns={"transcriptID": f"cat_{species}_transcriptID"})
    )
    checkPipe = collapse[collapse[f"cat_{species}_transcriptID"].str.contains(r"\|", na=False)]
    #print(checkPipe)
    print(f"rows in {species} cat: {len(collapse)}")
    print(f"# uniq cat in {species}:", collapse[f"cat_{species}_transcriptID"].nunique())
        
    # Sort transcript IDs inside the string
    collapse[f"cat_{species}_transcriptID"] = collapse[f"cat_{species}_transcriptID"].apply(
        lambda x: "|".join(sorted(x.split("|")))
    )
    anno[f"cat_{species}_transcriptID"] = anno[f"cat_{species}_transcriptID"].apply(
        lambda x: "|".join(sorted(x.split("|")))
    )

    ## merge collapse back to anno file
    anno_close = pd.merge(
        anno,
        collapse,
        on = f"{species}_jxnHash",
        how = 'outer',
        indicator = 'merge2')
    counts2 = anno_close["merge2"].value_counts()
    # Extract counts safely (0 if not present)
    right_only = counts2.get("right_only", 0)
    left_only = counts2.get("left_only", 0)
    both = counts2.get("both", 0)
    print(f"final merge for {species}: L_only: {left_only}, R_only: {right_only}, B: {both}, L+B: {left_only + both}")
        # L_only: 0, R_only: 0, B: 77023
    
    ##  keep if both
    anno_ready = anno_close[anno_close["merge2"].isin(["both"])].copy()
    ## confirm cat_x and cat_y are same
    mismatch = anno_ready[
        anno_ready[f"cat_{species}_transcriptID_x"] != anno_ready[f"cat_{species}_transcriptID_y"]
    ]
    # drop the merge indicators if you don’t need it
    anno_ready = anno_ready.drop(columns=["merge2"])
    print(anno_ready.columns)
    # rename flag column
    anno_ready = anno_ready.rename(columns={"flag_geneMismatch": "flag_GFFcomp_mismatch"})

    # check row counts before and after adding flag
    print(f"rows in {species} anno before: {len(anno)}")
    print(f"rows in {species} anno afer: {len(anno_ready)}")
    print(f"# uniq jxnHash in {species} after:", anno_ready[f"{species}_jxnHash"].nunique())
            
    ## save to csv
    anno_ready.to_csv(f"{proj}/submission/supplementary/fiveSpecies_{species}_full_annotation_w_dataFlag_gffCompMM.csv", index=False)    

