#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  1) replace incorrect geneIDs using the *_duplicated_geneID_fix.csv mapping
  2) update flag_jxnHash_in_fiveSpecies_full_anno = 1

"""

import pandas as pd

PROJ = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
DATA_DIR = f"{PROJ}/submission/supplementary_files/datafiles"
ANNO_DIR = f"{PROJ}/zenodo/fiveSpecies_supporting_files"


species_list = ['dmel6', 'dsim2', 'dyak2', 'dsan1', 'dser1']

for species in species_list:
    datafile = pd.read_csv(f"{DATA_DIR}/datafile_jxnHash_{species}_w_annoFlag_ERP_ESP_info_strCat_flagNovel_02amm.csv")

    if species in ['dmel6', 'dsim2']:
        # Fix 1: replace duplicated/missing geneIDs using mapping from Dr. Morse
        geneID_map_path = f"{ANNO_DIR}/{species}_duplicated_geneID_fix.csv"
        geneID_map = pd.read_csv(geneID_map_path, low_memory=False)
        
        id_mapping = dict(zip(geneID_map['lmm_missing_geneID'], geneID_map['geneID_NEW']))
        datafile['geneID'] = datafile['geneID'].replace(id_mapping)
        ## check one
        #check = datafile[datafile["geneID"] == "FBgn0032404"]
        #print((datafile["geneID"] == "FBgn0085193").any())
        
    # 2. check flag_jxnHash_in_fiveSpecies_full_anno for the jxnHashes
    anno = pd.read_csv(f"{PROJ}/zenodo/fiveSpecies_{species}_full_annotation.csv")
    keep_cols = [f"{species}_jxnHash"]
    anno = anno[keep_cols]
    
    datafile_anno = pd.merge(
        anno,
        datafile,
        on = f"{species}_jxnHash",
        how = 'outer',
        indicator = 'merge'
        )
    print(datafile_anno['merge'].value_counts(dropna=False).sort_index())
#            left_only       53314  anno only
#            right_only    2264256  datafile only
#            both            22672  in both
    
    ## flag_jxnHash_fiveSpecies_full_anno should be 1 if in both else 0
    datafile_anno["flag_jxnHash_fiveSpecies_full_anno"] = (
        datafile_anno["merge"] == "both"
    ).astype(int)

    #keep right and both, drop indicator col
    datafile_anno = datafile_anno[datafile_anno['merge'].isin(['right_only', 'both'])]
    
    datafile_anno.drop(columns="merge", inplace=True)
    print(datafile_anno["flag_jxnHash_fiveSpecies_full_anno"].value_counts())
    
    print(datafile_anno.columns)
    
    OUT_PATH = f"{DATA_DIR}/datafile_jxnHash_{species}_almost.csv"
    datafile_anno.to_csv(OUT_PATH, index=False)

