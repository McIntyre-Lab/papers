#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 00:18:44 2026

@author: nkeil
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# You need: pip install upsetplot
from upsetplot import UpSet, from_indicators

##############################################################################
# USER SETTINGS
##############################################################################

INPUT_DIR = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts"
OUTPUT_DIR = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

FILE_TEMPLATE = "metaanalysis_summary_{cell}_w_ES_diff_w_bysex.csv"
CELL_TYPES = ["CD4", "CD8"]

GENE_COL = "gene_id"
T1D_FLAG_COL = "flag_t1d_gene"
F_FLAG_COL = "flag_sig_heterogeneity_F"
M_FLAG_COL = "flag_sig_heterogeneity_M"

##############################################################################
# FUNCTIONS
##############################################################################

def load_and_check(path):
    df = pd.read_csv(path)

    required_cols = [GENE_COL, T1D_FLAG_COL, F_FLAG_COL, M_FLAG_COL]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    return df


def make_membership_count_table(df, cell, outdir):
    """
    Create a counts table for:
      - neither F nor M
      - F only
      - M only
      - both F and M
    Only for T1D genes.
    """

    df = df.copy()

    # Keep only T1D genes
    df = df[df[T1D_FLAG_COL] == 1].copy()

    # In case there are duplicate gene rows, keep one row per gene
    # If duplicates exist and any row has a 1, gene is treated as 1 for that flag
    df = (
        df.groupby(GENE_COL, as_index=False)[[F_FLAG_COL, M_FLAG_COL]]
          .max()
    )

    f = df[F_FLAG_COL].fillna(0).astype(int)
    m = df[M_FLAG_COL].fillna(0).astype(int)

    count_neither = ((f == 0) & (m == 0)).sum()
    count_f_only  = ((f == 1) & (m == 0)).sum()
    count_m_only  = ((f == 0) & (m == 1)).sum()
    count_both    = ((f == 1) & (m == 1)).sum()

    count_f_total = (f == 1).sum()
    count_m_total = (m == 1).sum()
    total_genes   = len(df)

    counts_df = pd.DataFrame({
        "cellType": [cell] * 7,
        "category": [
            "total_t1d_genes",
            "F_only",
            "M_only",
            "F_and_M",
            "neither",
            "F_total",
            "M_total"
        ],
        "count": [
            total_genes,
            count_f_only,
            count_m_only,
            count_both,
            count_neither,
            count_f_total,
            count_m_total
        ]
    })

    out_csv = os.path.join(outdir, f"upset_counts_{cell}_t1d_genes.csv")
    counts_df.to_csv(out_csv, index=False)

    print(f"[{cell}] Total T1D genes: {total_genes}")
    print(f"[{cell}] F only:  {count_f_only}")
    print(f"[{cell}] M only:  {count_m_only}")
    print(f"[{cell}] Both:    {count_both}")
    print(f"[{cell}] Neither: {count_neither}")
    print(f"[{cell}] Saved counts table: {out_csv}")

    return df, counts_df


def make_upset_plot(df, cell, outdir):
    """
    Make an upset plot from the two heterogeneity flags.
    """

    plot_df = df.copy()

    # Ensure boolean indicators
    plot_df[F_FLAG_COL] = plot_df[F_FLAG_COL].fillna(0).astype(int).astype(bool)
    plot_df[M_FLAG_COL] = plot_df[M_FLAG_COL].fillna(0).astype(int).astype(bool)

    # Create upsetplot input
    upset_data = from_indicators([F_FLAG_COL, M_FLAG_COL], plot_df)

    plt.figure(figsize=(8, 6))
    upset = UpSet(
        upset_data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality"
    )
    upset.plot()

    plt.suptitle(f"{cell}: T1D genes\nSignificant heterogeneity by sex", y=1.02)

    out_png = os.path.join(outdir, f"upset_{cell}_t1d_genes_heterogeneity_F_M.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[{cell}] Saved upset plot: {out_png}")


##############################################################################
# MAIN
##############################################################################

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for cell in CELL_TYPES:
        infile = os.path.join(INPUT_DIR, FILE_TEMPLATE.format(cell=cell))
        print(f"\nProcessing {cell}: {infile}")

        df = load_and_check(infile)
        df_t1d, counts_df = make_membership_count_table(df, cell, OUTPUT_DIR)
        make_upset_plot(df_t1d, cell, OUTPUT_DIR)


if __name__ == "__main__":
    main()