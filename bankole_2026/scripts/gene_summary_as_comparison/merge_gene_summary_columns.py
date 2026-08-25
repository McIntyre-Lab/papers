#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Outer-merge the gene-level UJC and ERP summary columns on geneID."""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Merge gene-level UJC and ERP summary columns into one table.")
    parser.add_argument("--ujc-summary", required=True)
    parser.add_argument("--erp-summary", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-l", "--log-file", required=True)
    args = parser.parse_args()

    ujc_df = pd.read_csv(args.ujc_summary, low_memory=False)
    erp_df = pd.read_csv(args.erp_summary, low_memory=False)

    gene_df = pd.merge(ujc_df, erp_df, on="geneID", how="outer", indicator="merge_check")

    with open(args.log_file, "w") as log:
        log.write(f"Merge gene summary columns: {args.species}\n")
        log.write(f"UJC columns: {ujc_df.shape[0]:,} genes, {ujc_df.shape[1]} columns\n")
        log.write(f"ERP columns: {erp_df.shape[0]:,} genes, {erp_df.shape[1]} columns\n")
        log.write("Merge ERP gene summary columns to UJC gene summary columns:\n")
        log.write(f"{gene_df['merge_check'].value_counts(dropna=False).sort_index()}\n")

    gene_df = gene_df.drop(columns=["merge_check"])
    gene_df = gene_df.fillna(0)
    gene_df.to_csv(args.output, index=False)

    with open(args.log_file, "a") as log:
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {gene_df.shape[0]:,} genes, {gene_df.shape[1]} columns\n")


if __name__ == "__main__":
    main()
