#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:06:18 2026

@author: ammorse


"""

import pandas as pd
import os
from scipy.stats import chi2
import numpy as np


proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

species_list = ["dmel6", "dsim2", "dser1"]
#species_list = ["dmel6"]
for species in species_list:

    # Import CSV
    infile = os.path.join(
        proj,
        f"{proj}/zenodo/datafiles/datafile_jxnHash_{species}.csv"
    )

    data = pd.read_csv(infile, low_memory = False)

    # Create exonERP
    data["exonERP"] = pd.NA


    ## for numExon_ERP thresholds of 2, 7, and 12
    set1 = data.copy()
    
    set1.loc[data["numExon_ERP"] >= 12, "exonERP"] = 3
    set1.loc[
        (set1["numExon_ERP"] >= 7) &
        (set1["numExon_ERP"] < 12),
        "exonERP"
    ] = 2
    set1.loc[
        (set1["numExon_ERP"] >= 2) &
        (set1["numExon_ERP"] < 7),
        "exonERP"
    ] = 1

    # Keep analyzable rows
    set1_2 = set1[set1["flag_analyzable"] == 1].copy()

    # tables exonERP * sexClass / out=crossTab;
    crossTab = pd.crosstab(
        set1_2["exonERP"],
        set1_2["sexClass"]
    ).reset_index()

    crossTab['bias'] = crossTab['F_bias'] + crossTab['M_bias']
    
    crossTab1_2 = crossTab.drop(columns=["F_bias", "M_bias", "no_ttest"])
    print(crossTab1_2)
    
    ## for numExon_ERP thresholds of 2, 4, and 7
    set2 = data.copy()
    
    set2.loc[set2["numExon_ERP"] >= 7, "exonERP"] = 3
    set2.loc[
        (set2["numExon_ERP"] >= 4) &
        (set2["numExon_ERP"] < 7),
        "exonERP"
    ] = 2
    set2.loc[
        (set2["numExon_ERP"] >= 2) &
        (set2["numExon_ERP"] < 4),
        "exonERP"
    ] = 1

    # Keep analyzable rows
    set2_2 = set2[set2["flag_analyzable"] == 1].copy()

    # tables exonERP * sexClass / out=crossTab;
    crossTab2 = pd.crosstab(
        set2_2["exonERP"],
        set2_2["sexClass"]
    ).reset_index()

    crossTab2['bias'] = crossTab2['F_bias'] + crossTab2['M_bias']
    
    crossTab2_2 = crossTab2.drop(columns=["F_bias", "M_bias", "no_ttest"])
    print(crossTab2_2)
    
    
    
    # Save
    crossTab1_2.to_csv(f"{proj}/submission/crossTab_exonERP_vs_sexClass_{species}_thr2_7_12.csv", index=False)
    
    crossTab2_2.to_csv(f"{proj}/submission/crossTab_exonERP_vs_sexClass_{species}_thr2_4_7.csv", index=False)
