#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:10:58 2026

@author: mgaran
"""
import numpy as np
import pandas as pd

PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
OUT = f"{PROJ}/zenodo/summary_files"
LOG_DIR = f"{PROJ}/Tables/log_files"

SPECIES_LIST = ["dmel", "dsim", "dser", "dyak", "dsan"]
SPECIES_10X = ["dyak", "dsan"]

G1 = "F"
G2 = "M"

def collapse_sexbias(df, bias_col, g1, g2):
    flags = (
        df.assign(
            has_B=df[bias_col].eq("B"),
            has_g1=df[bias_col].eq(g1),
            has_g2=df[bias_col].eq(g2),
            has_unbiased=df[bias_col].eq("unbiased"),
        )
        .groupby("genesetid", as_index=False)[["has_B", "has_g1", "has_g2", "has_unbiased"]]
        .max()
    )
    flags[bias_col] = np.select(
        [flags["has_B"], flags["has_g1"] & flags["has_g2"], flags["has_g1"], flags["has_g2"], flags["has_unbiased"]],
        ["B", "B", g1, g2, "unbiased"],
        default="not_evaluated",
    )
    return flags[["genesetid", bias_col]]


def collapse_sexbias_anno_status(df, anno_col):
    priority = ["anno_UJC_bias", "anno_ERP_bias", "novel_ERP_bias", "ISM_bias", "unbiased"]
    flags = (
        df.assign(**{f"has_{status}": df[anno_col].eq(status) for status in priority})
        .groupby("genesetid", as_index=False)[[f"has_{status}" for status in priority]]
        .max()
    )
    flags[anno_col] = np.select(
        [flags[f"has_{status}"] for status in priority],
        priority,
        default="not_evaluated",
    )
    return flags[["genesetid", anno_col]]


def main():
    for species in SPECIES_LIST:
        print(f"Processing {species}...")
        log_file = f"{LOG_DIR}/create_genesetid_summary_{species}.log"
        input_file = f"{OUT}/gene_level_summary_from_dataFile_jxnHash_w_genesetid_{species}.csv"
        output_file = f"{OUT}/genesetid_summary_{species}.csv"

        df = pd.read_csv(input_file, low_memory=False)

        with open(log_file, "w") as log:
            log.write(f"Create genesetid summary from gene summary: {species}\n\n")
            log.write(f"Input: {input_file}\n")
            log.write(f"Gene summary rows: {len(df):,}\n")
            log.write(f"Unique genes: {df['geneID'].nunique():,}\n")
            log.write(f"check genesetid is unique on geneID: {df.set_index('geneID')['genesetid'].is_unique}\n")
            n_no_gs = df["genesetid"].isna().sum()
            log.write(f"Genes with null genesetid: {n_no_gs:,}\n")
            if n_no_gs > 0:
                sample = df.loc[df["genesetid"].isna(), "geneID"].tolist()
                log.write(f"  geneIDs: {sample}\n")
            log.write("\n")

        gs_df = df.dropna(subset=["genesetid"]).copy()
        gs_df["genesetid"] = gs_df["genesetid"].astype(int)

        agg_cols = {
            "num_geneID": ("geneID", "nunique"),
            "sumReads": ("sumReads", "sum"),
            f"sumReads_{G1}": (f"sumReads_{G1}", "sum"),
            f"sumReads_{G2}": (f"sumReads_{G2}", "sum"),
        }
        for col in ["num_components", "num_components_simple"]:
            if col in gs_df.columns:
                agg_cols[col] = (col, "sum")

        summary = gs_df.groupby("genesetid").agg(**agg_cols).reset_index()

        if species in SPECIES_10X:
            bias_col = "sexBias_10x"
            anno_col = "sexBias_10x_anno_status"
            ujcan_col = "sexBias_10x_UJCanno"
        else:
            bias_col = "sexBias"
            anno_col = "sexBias_anno_status"
            ujcan_col = "sexBias_UJCanno"

        for col in [bias_col, anno_col, ujcan_col]:
            if col not in gs_df.columns:
                with open(log_file, "a") as log:
                    log.write(f"WARNING: column '{col}' not found in gene summary; skipping\n")
                continue

            collapsed = collapse_sexbias_anno_status(gs_df, col) if col == anno_col else collapse_sexbias(gs_df, col, G1, G2)
            summary = pd.merge(summary, collapsed, on="genesetid", how="outer", indicator="merge_check")
            with open(log_file, "a") as log:
                log.write(f"Merge {col} to summary:\n")
                log.write(f"{summary['merge_check'].value_counts(dropna=False).sort_index()}\n\n")
            summary = summary[summary["merge_check"] != "right_only"].drop(columns=["merge_check"])

        summary.to_csv(output_file, index=False)

        with open(log_file, "a") as log:
            log.write(f"Output: {output_file}\n")
            log.write(f"Final: {len(summary):,} genesetids, {summary.shape[1]} columns\n")

        print(f"  {len(summary):,} genesetids -> {output_file}")

if __name__ == "__main__":
    main()