#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 13:47:14 2025

@author: nkeil
"""

import pandas as pd
import os 
import numpy as np

ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/datafiles"
outd= "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis"

# Set up tag directory
species_tags = {
    "mel": "dmel6",
    "sim": "dsim2",
    "yak": "dyak2",
    "san": "dsan1",
    "ser": "dser1"
}

#Set up columns that we want to keep
base_keep_cols = [
    "geneID",
    "ERP",
    "{tag}_jxnHash",
    "structural_category",
    "flag_ERPanno",
    "flag_putative_novel_jxnHash",
]

#Structural categpry abbreviations
abbr_mapping = {
    "full-splice_match": "FSM",
    "incomplete-splice_match": "ISM",
    "novel_in_catalog": "NIC",
    "novel_not_in_catalog": "NNC",
    "antisense": "AS",
    "fusion": "FUS",
    "genic": "GENIC",
    "genic_intron": "GI",
    "intergenic": "INTER",
    "ambiguous": "AMBIG"
}

#Set up cnts prefix
cnts_prefix = "rawCnts"  # everything starting with this will be kept

#Function to make list of columns to keep
def make_usecols(keep_cols, prefix=cnts_prefix):
    """Return a callable for pd.read_csv(usecols=...) that keeps selected columns."""
    def _usecols(colname: str) -> bool:
        return colname in keep_cols or colname.startswith(cnts_prefix)
    return _usecols

for species, tag in species_tags.items():
    # expand the {tag} placeholder
    keep_cols = [c.format(tag=tag) if "{tag}" in c else c for c in base_keep_cols]

    datafile = f"datafile_jxnHash_{tag}_w_annoFlag_ERP_ESP_info_strCat_flagNovel.csv"
    datafile_path = os.path.join(ind, datafile)

    # load ONLY required columns
    data_df = pd.read_csv(
        datafile_path,
        usecols=make_usecols(keep_cols, cnts_prefix),
        low_memory=False,  # safer mixed-type parsing
    )
    
    #Create structural category plus 
    # Map to abbreviations first
    data_df["structural_category_abbr"] = data_df["structural_category"].map(abbr_mapping)
    
    # Apply rules for the new column
    data_df["structural_category_plus"] = np.where(
        data_df["structural_category_abbr"] == "NNC",
        np.where(data_df["flag_ERPanno"] == 1,
                 "NNC_annoERP",
                 "NNC_no_annoERP"),
        data_df["structural_category_abbr"]
    )
    
    data_df.drop(columns=["structural_category_abbr"], inplace=True)
    
    data_df.to_csv(os.path.join(outd,f"datafile_jxnHash_{tag}_w_struct_category_plus.csv"), index=False)


    
    