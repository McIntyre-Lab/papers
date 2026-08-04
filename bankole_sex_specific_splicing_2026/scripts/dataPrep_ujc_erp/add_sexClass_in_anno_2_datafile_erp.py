#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 11:09:36 2025

@author: ammorse
"""

import pandas as pd
import numpy as np

proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

# jxnHash sexClass for mel sim and ser:
species_dct = {"dmel6":"dmel", "dsim2":"dsim", "dser1":"dser"}
#species_dct = {"dmel6":"dmel"}

for species,species_data in species_dct.items():
    print(f"processing species: {species}")    
    # import datafile                            
    datafile = pd.read_csv(f"{proj}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag_flagNovel.csv", low_memory=False)
                           
    print(datafile.columns)
    cols_with_pval = [c for c in datafile.columns if "pval" in c.lower()]
    print(cols_with_pval)

    # where flag_in_full_anno = 1 create sexClass
    tmp = np.select(
        [
            (datafile.flag_sex_limited_F_rc50 == 1),
            (datafile.flag_sex_limited_M_rc50 == 1),
            (datafile.flag_analyzable == 0),
            (datafile.flag_sex_bias_F == 1),
            (datafile.flag_sex_bias_M == 1),
            (datafile.pval_ttest_Equal > 0.05),
        ],
        [
            "F_limited",
            "M_limited",
            "unanalyzed",
            "F_bias",
            "M_bias",
            "unbiased",
        ],
        default="no_ttest"
    )
    
    mask = datafile.flag_in_full_anno == 1
    datafile.loc[mask, "sexClass_in_anno"] = tmp[mask]
    
    datafile.loc[datafile.flag_in_full_anno == 0, "sexClass_in_anno"] = np.nan
    
    # check with crosstab
    freq = pd.crosstab(datafile['sexClass_in_anno'], datafile['flag_in_full_anno'])
    #print(freq)
    ones = freq[1.0].sum()
    print(ones)  ## 32487    
    print(datafile.flag_in_full_anno.value_counts()) ## 32487
    
    # put sexClass_in_anno after sexClass
    cols = list(datafile.columns)
    if "sexClass_in_anno" in cols:
        cols.remove("sexClass_in_anno")
    # find index of 'sexClass' and insert new column right after
    idx = cols.index("sexClass")
    cols.insert(idx + 1, "sexClass_in_anno")
    # reorder the DataFrame
    datafile = datafile[cols]

    # save to csv
    datafile.to_csv(f"{proj}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag_flagNovel_02amm.csv", index=False)    
    
    
# jxnHash sexClass for yak and san (no pvals) where flag_anno = 1
species_dct = {"dyak2":"dyak", "dsan1":"dsan"}
#species_dct = {"dyak2":"dyak"}

for species,species_data in species_dct.items():
    print(f"processing species: {species}")       
    # import datafile                            
    datafile = pd.read_csv(f"{proj}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag_flagNovel.csv")
    print(datafile.columns)

    # where flag_jxnHash_in_fiveSpecies_full_anno = 1 create sexClass
    tmp = np.select(
        [
            (datafile.flag_sex_limited_F_rc50 == 1),
            (datafile.flag_sex_limited_M_rc50 == 1),
            (datafile.flag_analyzable == 0),
            (datafile.flag_sex_bias_F == 1),
            (datafile.flag_sex_bias_M == 1),
        ],
        [
            "F_limited",
            "M_limited",
            "unanalyzed",
            "F_bias",
            "M_bias",
        ],
        default="unanalyzed"
    )
    
    mask = datafile.flag_in_full_anno == 1
    datafile.loc[mask, "sexClass_in_anno"] = tmp[mask]

    datafile.loc[datafile.flag_in_full_anno == 0, "sexClass_in_anno"] = np.nan
    
    # check with crosstab
    freq = pd.crosstab(datafile['sexClass_in_anno'], datafile['flag_in_full_anno'])
    #print(freq)
    ones = freq[1.0].sum()
    print(ones)  ## 17937
    
    print(datafile.flag_in_full_anno.value_counts()) ## 17937

    # put sexClass_in_anno after sexClass
    cols = list(datafile.columns)
    if "sexClass_in_anno" in cols:
        cols.remove("sexClass_in_anno")
    # find index of 'sexClass' and insert new column right after
    idx = cols.index("sexClass")
    cols.insert(idx + 1, "sexClass_in_anno")
    # reorder the DataFrame
    datafile = datafile[cols]
    
    # save to csv
    datafile.to_csv(f"{proj}/submission/supplementary/datafiles/datafile_erp_{species}_w_annoFlag_flagNovel_02amm.csv", index=False)    
    