#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge genesetID into gene summaries and place it as the second column.

Input/output file stems are given as arguments so the same script serves the
full and min10reads runs; per-species path is f"{summary_dir}/{prefix}_{species}.csv".
"""

import argparse
import pandas as pd

PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
NODES_FILE = f"{PROJ}/submission/supplementary_files/fiveSpecies_annotations/link_files/nodes_with_geneset.csv"
LOG = f"{PROJ}/Tables/log_files/merge_genesetID_to_geneSymbol_gene_summaries.log"

SPECIES_LIST = ["dmel", "dsim", "dser", "dyak", "dsan"]


def main():
    parser = argparse.ArgumentParser(description="Merge genesetID into gene summaries.")
    parser.add_argument("--summary-dir", default=f"{PROJ}/zenodo/summary_files")
    parser.add_argument("--in-prefix", default="gene_level_summary_w_geneSymbol")
    parser.add_argument("--out-prefix", default="gene_summary")
    args = parser.parse_args()

    OUT = args.summary_dir

    nodes_df = pd.read_csv(NODES_FILE, dtype=str)
    gene_nodes = nodes_df[~nodes_df["nodeid"].astype(str).str.startswith("c")].copy()
    gene_nodes = gene_nodes.rename(columns={"nodeid": "geneID"})
    gene_nodes = gene_nodes[["geneID", "genesetid"]].drop_duplicates()
    gene_nodes["genesetid"] = gene_nodes["genesetid"].astype(int)

    with open(LOG, "w") as log:
        log.write("Merge genesetID to gene summaries after gene symbol merge\n\n")
        log.write(f"Gene nodes available for merge: {len(gene_nodes):,}\n")
        log.write(f"Check gene nodes unique on geneID: {gene_nodes['geneID'].is_unique}\n\n")

    for species in SPECIES_LIST:
        input_file = f"{OUT}/{args.in_prefix}_{species}.csv"
        output_file = f"{OUT}/{args.out_prefix}_{species}.csv"

        gene_df = pd.read_csv(input_file, low_memory=False)
        if "genesetid" in gene_df.columns:
            gene_df = gene_df.drop(columns=["genesetid"])

        merged = pd.merge(gene_df, gene_nodes, on="geneID", how="left", indicator="merge_check")

        with open(LOG, "a") as log:
            log.write(f"{species}:\n")
            log.write(f"  Gene summary rows: {len(gene_df):,}\n")
            log.write("  Merge results:\n")
            log.write(f"{merged['merge_check'].value_counts(dropna=False).sort_index()}\n")
            log.write(f"  Rows with genesetid: {merged['genesetid'].notna().sum():,}\n")
            log.write(f"  Rows without genesetid: {merged['genesetid'].isna().sum():,}\n\n")

        merged = merged.drop(columns=["merge_check"])
        cols = merged.columns.tolist()
        cols.remove("genesetid")
        cols.insert(1, "genesetid")
        merged = merged[cols]
        merged.to_csv(output_file, index=False)
        print(f"{species}: saved {len(merged):,} rows to {output_file}")

    print(f"Log: {LOG}")


if __name__ == "__main__":
    main()
