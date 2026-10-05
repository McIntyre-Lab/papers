#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 13:32:08 2026

@author: ammorse
"""


import pandas as pd
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
import numpy as np


proj = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

species_list = ["dmel6", "dsim2", "dser1"]
#species_list = ["dmel6"]

results = []

for species in species_list:

    # Import CSV
    data = pd.read_csv(f"{proj}/zenodo/datafiles/datafile_jxnHash_{species}.csv", low_memory = False)
    
    # Keep analyzable rows
    data2 = data[data["flag_analyzable"] == 1].copy()
    
    ## compare exons in biased ujc vs in unbiased ujc
    data2["bias_status"] = np.where(
    data2["sexClass"].isin(["F_bias", "M_bias"]),
    "bias",
    np.where(
        data2["sexClass"] == "unbiased",
        "unbiased",
        np.nan
        )
    )
    
    ch = data2["bias_status"].value_counts()
    
    biased = data2.loc[
        data2["bias_status"] == "bias",
        "numExon_ERP"
    ].dropna()
    
    unbiased = data2.loc[
        data2["bias_status"] == "unbiased",
        "numExon_ERP"
    ].dropna()
    
    tstat, pvalue_t = ttest_ind(
        biased,
        unbiased,
        equal_var=False
    )

    print("Biased UJCs:", len(biased))
    print("Unbiased UJCs:", len(unbiased))
    print("Mean exons - biased:", biased.mean())
    print("Mean exons - unbiased:", unbiased.mean())
    print("t-statistic:", tstat)
    print("p-value:", pvalue_t)

#    Biased UJCs: 10719
#    Unbiased UJCs: 49397
#    Mean exons - biased: 3.014087134993936
#    Mean exons - unbiased: 3.0794785108407394
#    t-statistic: -3.3790576230024993
#    p-value: 0.0007289088339415467

    u, pvalue_mw = mannwhitneyu(
        biased,
        unbiased,
        alternative="two-sided"
    )
    
    print("Mann-Whitney U:", u)
    print("p-value:", pvalue_mw)

#    Mann-Whitney U: 271264308.5
#    p-value: 4.364623390907274e-05

    print("Biased:")
    print("n =", len(biased))
    print("mean =", biased.mean())
    print("median =", biased.median())
    print("SD =", biased.std())
    
    print("\nUnbiased:")
    print("n =", len(unbiased))
    print("mean =", unbiased.mean())
    print("median =", unbiased.median())
    print("SD =", unbiased.std())

    # Biased:
    # n = 10719
    # mean = 3.014087134993936
    # median = 3.0
    # SD = 1.7566251008601967
    
    # Unbiased:
    # n = 49397
    # mean = 3.0794785108407394
    # median = 3.0
    # SD = 2.0685571800941958

    # Save results
    results.append({
        "species": species,
        "biased_n": len(biased),
        "unbiased_n": len(unbiased),
        "biased_mean": biased.mean(),
        "unbiased_mean": unbiased.mean(),
        "biased_median": biased.median(),
        "unbiased_median": unbiased.median(),
        "biased_SD": biased.std(),
        "unbiased_SD": unbiased.std(),
        "t_stat": tstat,
        "t_pvalue": pvalue_t,
        "MW_U": u,
        "MW_pvalue": pvalue_mw
    })
    
# Results table
results_df = pd.DataFrame(results)

results_df.to_csv(f"{proj}/submission/ttest_numExon_ERP_vs_bias.csv", index=False)
    
print("\n\nSUMMARY")
print(results_df)    
    
    