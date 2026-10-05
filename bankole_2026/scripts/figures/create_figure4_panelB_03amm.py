#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

BASE_DIR = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
FIGURES_DIR = f"{BASE_DIR}/Figures_amm"
SUMMARY_DIR = f"{BASE_DIR}/zenodo/summary_files"
DUP_DIR = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/duplication_vs_splicing/analysis_adp"

dark = {"dmel": "#966729", "dsim": "#3F78C1", "dsan": "#28827A", "dyak": "#717273", "dser": "#825CA6"}

rng = np.random.default_rng(seed=42)

panel_b = [
    ("dmel", f"{DUP_DIR}/dmel_singleton_vs_dup.csv", {"dmel": "geneID"}, {}, "geneID", "num_FBGNs", 0.15, "Melanogaster"),
    ("dsim", f"{DUP_DIR}/dsim_5species_Fbgn_2_Hahn_GLEANR.csv", {}, {"geneID": "dsim_FBgn"}, "dsim_FBgn", "num_GLEANRIDs_sim", 0.001, "Simulans"),
]

row_specs = [
    ("w", "Genes with sex bias"),
    ("wo", "Genes without statistical evidence for sex bias (unbiased)"),
]

fig, axes = plt.subplots(
    2, 4, figsize=(18, 12), sharey=True,
    gridspec_kw={"width_ratios": [2, 1, 2, 1], "wspace": 0.08}
)

for col, (sp, dup_path, dup_rename, data_rename, merge_on, xcol, jit, name) in enumerate(panel_b):
    dup = pd.read_csv(dup_path).rename(columns=dup_rename)
    data = pd.read_csv(f"{SUMMARY_DIR}/gene_summary_{sp}.csv", low_memory=False).rename(columns=data_rename)
    #print(data.columns[91])
    #print(data.iloc[:, 91].unique())


    data["num_nonISM_ERPp"] = (
        data["num_ERPp_w_anno_ujc"]
        + data["num_ERPp_w_anno_erp"]
        + data["num_ERPp_w_novel_erp"]
    )
    
    # nan check - should be all 0
    print(data[["num_ERPp_w_anno_ujc","num_ERPp_w_anno_erp","num_ERPp_w_novel_erp"]].isna().sum())
    
    merged = pd.merge(dup, data, on=merge_on, how="inner")
    x = merged[xcol].to_numpy(dtype=float)
    y = merged["num_nonISM_ERPp"].to_numpy(dtype=float)

    #merged["sexBias"].value_counts(dropna=False)
    
    #mask = merged["sexBias"].isin(["M", "F", "B"]).to_numpy()
    mask_biased = merged["sexBias"].isin(["M", "F", "B"]).to_numpy()
    mask_unbiased = merged["sexBias"].isin(["unbiased"]).to_numpy()
    
    left_idx = col * 2
    right_idx = col * 2 + 1

    for r, (sub, tag) in enumerate(row_specs):
        if sub == "w":
            m = mask_biased
        else:
            m = mask_unbiased
            
        #m = mask if sub == "w" else ~mask   # ~mask is not_evaluated and unbiased
        print(f"{sp}")
        print(merged["sexBias"].unique())
        print(f"  w  (M/F/B):  {mask_biased.sum()}")
        print(f"  wo (unbiased):  {(mask_unbiased).sum()}")
        

        ax_l = axes[r, left_idx]
        ax_r = axes[r, right_idx]

        x_val = x[m]
        y_val = y[m]
       
        m_l = x_val <= 50   # broken x axis
        m_r = x_val >= 150  # broken x axis

        ax_l.scatter(
            x_val[m_l] + rng.uniform(-jit, jit, m_l.sum()), y_val[m_l],
            s=10, alpha=0.4, c=dark[sp], edgecolors="none", linewidths=0
        )
        ax_r.scatter(
            x_val[m_r] + rng.uniform(-jit, jit, m_r.sum()), y_val[m_r],
            s=10, alpha=0.4, c=dark[sp], edgecolors="none", linewidths=0
        )

        print(f"{sp} {sub}: max num_nonISM_ERPp = {y_val.max()}")
        #dmel w: max num_nonISM_ERPp = 373.0
        #dmel wo: max num_nonISM_ERPp = 239.0
        #dsim w: max num_nonISM_ERPp = 539.0    
        #dsim wo: max num_nonISM_ERPp = 308.0
        
        ax_l.set_xlim(0, 50)
        ax_r.set_xlim(150, 250)
        ax_l.set_ylim(0, 560)   ## amm increased so no clipping
        ax_r.set_ylim(0, 560)   ## amm increased so no clipping

        ax_l.spines["right"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        ax_l.yaxis.tick_left()
        ax_r.yaxis.tick_right()
        ax_r.tick_params(labelleft=False)

        if r == 1:
            ax_l.set_xlabel("Number of genes in gene family")
            ax_r.set_xlabel("")
        if left_idx == 0 and r == 0:
            ax_l.set_ylabel("number of ERP_plus without ISM reads")

        ax_l.set_title(f"{name} — {tag}", loc="left", fontsize=11, fontweight="bold")

        d = 0.015
        kwargs = dict(transform=ax_l.transAxes, color="k", clip_on=False, linewidth=1)
        ax_l.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax_l.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        kwargs.update(transform=ax_r.transAxes)
        ax_r.plot((-d, +d), (-d, +d), **kwargs)
        ax_r.plot((-d, +d), (1 - d, 1 + d), **kwargs)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/figure4_panelB_02mdg.pdf", format="pdf", bbox_inches="tight", dpi=600)
plt.close()
