#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create an ERP_plus-level propReads table.

Rolls jxnHash reads to ERP_plus (all/F/M), records per-ERP_plus counts of
jxnHash by annotation partition (anno_ujc / w_anno_erp / w_novel_erp / ism),
computes gene-normalized read proportions, and marks the most expressed
ERP_plus per gene.
"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Create ERP_plus-level propReads table.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("-p", "--count-prefix", required=True)
    parser.add_argument("--jxn-hash-col", default="jxnHash")
    parser.add_argument("--ujc-annotation-flag", default="flag_jxnHash_in_fiveSpecies_full_anno")
    parser.add_argument("--erp-annotation-flag", default="flag_ERPanno")
    args = parser.parse_args()

    g1 = args.group1
    g2 = args.group2
    anno_flag = args.ujc_annotation_flag
    erp_flag = args.erp_annotation_flag

    df = pd.read_csv(args.input, low_memory=False)
    df = df.rename(columns={args.jxn_hash_col: "jxnHash"})

    total_count_columns = [column for column in df.columns if column.startswith(f"{args.count_prefix}_")]
    group1_count_columns = [column for column in total_count_columns if f"_{g1}_" in column]
    group2_count_columns = [column for column in total_count_columns if f"_{g2}_" in column]
    df[total_count_columns] = df[total_count_columns].apply(pd.to_numeric, errors="coerce").fillna(0)

    df["sumReads_T"] = df[total_count_columns].sum(axis=1)
    df[f"sumReads_{g1}"] = df[group1_count_columns].sum(axis=1)
    df[f"sumReads_{g2}"] = df[group2_count_columns].sum(axis=1)

    df["sumReads_anno_ujc"] = (df[anno_flag] == 1) * df["sumReads_T"]
    df["sumReads_anno_erp"] = (df[erp_flag] == 1) * df["sumReads_T"]
    df["sumReads_novel_erp"] = (df[erp_flag] == 0) * df["sumReads_T"]
    df["sumReads_ism"] = (df["structural_category"] == "incomplete-splice_match") * df["sumReads_T"]

    erp_df = df.groupby(["geneID", "ERP_plus"], sort=False).agg(
        sumReads_T=("sumReads_T", "sum"),
        **{f"sumReads_{g1}": (f"sumReads_{g1}", "sum")},
        **{f"sumReads_{g2}": (f"sumReads_{g2}", "sum")},
        sumReads_anno_ujc=("sumReads_anno_ujc", "sum"),
        sumReads_anno_erp=("sumReads_anno_erp", "sum"),
        sumReads_novel_erp=("sumReads_novel_erp", "sum"),
        sumReads_ism=("sumReads_ism", "sum"),
        ERP=("ERP", "first"),
    ).reset_index()

    erp_df["numExon_GM"] = erp_df["ERP"].astype(str).str.len() - 2

    # For each ERP_plus, calculate read proportions
    for group in ["T", g1, g2]:
        reads = erp_df[f"sumReads_{group}"]
        gene_total = reads.groupby(erp_df["geneID"]).transform("sum")
        erp_df[f"{group}_read_proportion"] = (reads / gene_total).where(gene_total > 0, 0)

    # determine the most expressed ERP_plus per gene
    for group in ["T", g1, g2]:
        reads = erp_df[f"sumReads_{group}"]

        # Flag most expressed ERP_plus
        gene_max = reads.groupby(erp_df["geneID"]).transform("max")
        erp_df[f"flag_mstExpr_ERPp_{group}"] = ((reads == gene_max) & (gene_max > 0)).astype(int)

    erp_df.to_csv(args.output, index=False)

    print(f"Species/dataset: {args.species}", flush=True)
    print(f"Groups: {g1}, {g2}", flush=True)
    print(f"Genes: {erp_df['geneID'].nunique():,}", flush=True)
    print(f"ERP_plus rows: {erp_df.shape[0]:,}", flush=True)
    print(f"Output saved to: {args.output}", flush=True)

if __name__ == "__main__":
    main()