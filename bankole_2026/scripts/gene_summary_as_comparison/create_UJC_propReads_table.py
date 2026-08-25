#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a UJC-level propReads table.

Sums reads per jxnHash (all/F/M samples), rolls gene totals, computes read
proportions (all-UJC and annotated-UJC denominators), and marks the most
expressed jxnHash per gene across all UJC and across annotated UJC only.
Gene-level normalization is vectorized with groupby transforms.
"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Create UJC-level propReads table.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g", "--genome", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("-p", "--count-prefix", required=True)
    parser.add_argument("--ujc-annotation-flag", default="flag_jxnHash_in_fiveSpecies_full_anno")
    parser.add_argument("--erp-annotation-flag", default="flag_ERPanno")
    args = parser.parse_args()

    group1 = args.group1
    group2 = args.group2
    anno_flag = args.ujc_annotation_flag

    df = pd.read_csv(args.input, low_memory=False)
    df = df.rename(columns={f"{args.genome}_jxnHash": "jxnHash"})

    total_count_columns = [column for column in df.columns if column.startswith(f"{args.count_prefix}_")]
    group1_count_columns = [column for column in total_count_columns if f"_{group1}_" in column]
    group2_count_columns = [column for column in total_count_columns if f"_{group2}_" in column]

    df[total_count_columns] = df[total_count_columns].apply(pd.to_numeric, errors="coerce").fillna(0)

    df["sumReads_T"] = df[total_count_columns].sum(axis=1)
    df[f"sumReads_{group1}"] = df[group1_count_columns].sum(axis=1)
    df[f"sumReads_{group2}"] = df[group2_count_columns].sum(axis=1)
    df["sumReads_monoExon"] = df["sumReads_T"] * (df["numExon_ERP"] == 1)

    df["flag_is_incomplete-splice_match"] = (df["structural_category"] == "incomplete-splice_match").astype(int)

    # For each UJC, calculate read proportions
    for group_name in ["T", group1, group2]:
        reads = df[f"sumReads_{group_name}"]
        gene_total = reads.groupby(df["geneID"]).transform("sum")
        df[f"{group_name}_read_proportion"] = (reads / gene_total).where(gene_total > 0, 0)

        anno_reads = reads.where(df[anno_flag] == 1, 0)
        gene_total_anno = anno_reads.groupby(df["geneID"]).transform("sum")
        df[f"{group_name}_read_proportion_anno"] = (reads / gene_total_anno).where(
            df[anno_flag] == 1 & (gene_total_anno > 0), 0
        )

    # determine the most expressed UJC and annotated UJC per gene
    for group_name in ["T", group1, group2]:
        reads = df[f"sumReads_{group_name}"]
        anno_reads = reads.where(df[anno_flag] == 1)

        # Flag most expressed UJC
        gene_max = reads.groupby(df["geneID"]).transform("max")
        df[f"flag_mstExpr_UJC_{group_name}"] = ((reads == gene_max) & (gene_max > 0)).astype(int)
    
        # Flag most expressed annotated UJC
        gene_max_anno = anno_reads.groupby(df["geneID"]).transform("max")
        df[f"flag_mstExpr_annoUJC_{group_name}"] = ((anno_reads == gene_max_anno) & (gene_max_anno > 0)).astype(int)

    df = df.drop(columns=["structural_category"] + total_count_columns)
    df.to_csv(args.output, index=False)

    print(f"Species/dataset: {args.species}", flush=True)
    print(f"Groups: {group1}, {group2}", flush=True)
    print(f"Count columns: {len(total_count_columns):,}", flush=True)
    print(f"Genes: {df['geneID'].nunique():,}", flush=True)
    print(f"UJC rows: {df.shape[0]:,}", flush=True)
    print(f"Output saved to: {args.output}", flush=True)

if __name__ == "__main__":
    main()