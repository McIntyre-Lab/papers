#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

species_list = ["dmel", "dsim", "dsan", "dyak", "dser"]
genome_version = {"dmel": "dmel6", "dsim": "dsim2", "dsan": "dsan1", "dyak": "dyak2", "dser": "dser1"}

geneSymbol_path = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase650/dmel_annotation/fbgn_annotation_ID.csv"

PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
DATA = f"{PROJ}/submission/supplementary_files/datafiles"
log_dir = f"{PROJ}/Tables/log_files"


def main():
    parser = argparse.ArgumentParser(description="Merge dmel6 ortholog geneID and gene symbol into gene summaries.")
    parser.add_argument("--summary-dir", default=f"{PROJ}/zenodo/summary_files")
    parser.add_argument("--in-prefix", default="gene_level_summary_w_flags")
    parser.add_argument("--out-prefix", default="gene_level_summary_w_geneSymbol")
    args = parser.parse_args()

    OUT = args.summary_dir

    dmel_symbol_df = pd.read_csv(geneSymbol_path, usecols=["symbol", "primary_FBgn"], low_memory=False)
    dmel_symbol_df = dmel_symbol_df.rename(columns={"primary_FBgn": "dmel6_geneID", "symbol": "dmel_symbol"})

    for species in species_list:
        print(f"{'='*70}\nProcessing {species}\n{'='*70}", flush=True)
        genome = genome_version[species]
        input_df = f"{OUT}/{args.in_prefix}_{species}.csv"
        output_df = f"{OUT}/{args.out_prefix}_{species}.csv"
        log_file = f"{log_dir}/merging_dmel_orthologs_{args.out_prefix}_{species}.log"

        if species == "dmel":
            gene_df = pd.read_csv(input_df, low_memory=False)
            gene_df = gene_df.rename(columns={"geneID": "dmel6_geneID"})
            merged_df = pd.merge(gene_df, dmel_symbol_df, on="dmel6_geneID", how="outer", indicator="merge_check")
            with open(log_file, "w") as log:
                log.write("Merge dmel symbol to dmel geneID:\n")
                log.write(f"{merged_df['merge_check'].value_counts(dropna=False).sort_index()}\n")
            merged_df = merged_df[merged_df["merge_check"].isin(["left_only", "both"])].drop(columns=["merge_check"])
            merged_df = merged_df.rename(columns={"dmel6_geneID": "geneID"})
            merged_df.to_csv(output_df, index=False)
            print(f"Saved {len(merged_df):,} genes to {output_df}", flush=True)
            continue

        jxnHash_table_path = f"{DATA}/datafile_jxnHash_{genome}_w_annoFlag_ERP_ESP_info_strCat_flagNovel_02amm.csv"
        jxnHash_df = pd.read_csv(jxnHash_table_path, usecols=["geneID", "dmel6_geneID"], low_memory=False)
        ortholog_df = pd.merge(jxnHash_df, dmel_symbol_df, on="dmel6_geneID", how="outer", indicator="merge_check").drop_duplicates()
        with open(log_file, "w") as log:
            log.write(f"Merge dmel symbol to dmel6 ortholog column ({species}):\n")
            log.write(f"{ortholog_df['merge_check'].value_counts(dropna=False).sort_index()}\n")
        ortholog_df = ortholog_df[ortholog_df["merge_check"].isin(["left_only", "both"])].drop(columns=["merge_check"])

        gene_df = pd.read_csv(input_df, low_memory=False)
        merged_df = pd.merge(gene_df, ortholog_df, on="geneID", how="outer", indicator="merge_check")
        final_df = merged_df.groupby("geneID").agg({
            **{col: "first" for col in gene_df.columns if col != "geneID"},
            "dmel6_geneID": lambda x: "|".join(x.dropna().astype(str)),
            "dmel_symbol": lambda x: "|".join(x.dropna().astype(str)) if x.dropna().any() else "",
            "merge_check": "first",
        }).reset_index()
        with open(log_file, "a") as log:
            log.write("Merge ortholog list to gene summary:\n")
            log.write(f"{merged_df['merge_check'].value_counts(dropna=False).sort_index()}\n")
        final_df = final_df[final_df["merge_check"].isin(["left_only", "both"])].drop(columns=["merge_check"])
        final_df.to_csv(output_df, index=False)
        print(f"Saved {len(final_df):,} genes to {output_df}", flush=True)

    print("All species processing complete!", flush=True)


if __name__ == "__main__":
    main()
