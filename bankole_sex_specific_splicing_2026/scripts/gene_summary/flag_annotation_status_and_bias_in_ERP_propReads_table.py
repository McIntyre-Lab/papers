#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Flag annotation presence and bias in an ERP propReads table.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("-l", "--log-file", required=True)
    parser.add_argument("--erp-data-file", required=True)
    parser.add_argument("--group-class-col", default="sexClass")
    parser.add_argument("--analyzable-col", default="flag_analyzable")
    args = parser.parse_args()

    g1 = args.group1
    g2 = args.group2

    erp_df = pd.read_csv(args.input, low_memory=False)

    # Mutually exclusive annotation status assigned hierarchically: anno_UJC > ISM > anno_ERP > novel_ERP),
    erp_df["annotation_status"] = np.select(
        [erp_df["sumReads_anno_ujc"] >= 1, erp_df["sumReads_ism"] >= 1, erp_df["sumReads_anno_erp"] >= 0],
        ["anno_UJC", "ISM", "anno_ERP"],
        default="novel_ERP",
    )
    erp_df["flag_erpp_w_anno_ujc"] = (erp_df["annotation_status"] == "anno_UJC").astype(int)
    erp_df["flag_erpp_w_anno_erp"] = (erp_df["annotation_status"] == "anno_ERP").astype(int)
    erp_df["flag_erpp_w_novel_erp"] = (erp_df["annotation_status"] == "novel_ERP").astype(int)
    erp_df["flag_erpp_w_ism"] = (erp_df["annotation_status"] == "ISM").astype(int)

    stat_df = pd.read_csv(args.erp_data_file, low_memory=False)
    stat_df = stat_df[["geneID", "ERP_plus", args.group_class_col, args.analyzable_col]].copy()
    stat_df = stat_df.rename(columns={args.group_class_col: "sexClass", args.analyzable_col: "flag_analyzable"})
    stat_df = stat_df.drop_duplicates(subset=["geneID", "ERP_plus"])

    erp_df = pd.merge(erp_df, stat_df, on=["geneID", "ERP_plus"], how="left", indicator="merge_check")
    with open(args.log_file, "w") as log:
        log.write(f"Flag ERP propReads table: {args.species}\n")
        log.write(f"Input: {args.input}\n")
        log.write(f"ERP datafile: {args.erp_data_file}\n")
        log.write("Merge ERP statistical datafile to ERP propReads table:\n")
        log.write(f"{erp_df['merge_check'].value_counts(dropna=False).sort_index()}\n")
    erp_df = erp_df.drop(columns=["merge_check"])

    erp_df["sexClass"] = erp_df["sexClass"].fillna("")
    erp_df["flag_analyzable"] = erp_df["flag_analyzable"].fillna(0).astype(int)

    analyzable = erp_df["flag_analyzable"] == 1
    erp_df[f"flag_erpp_bias_{g1}"] = (erp_df["sexClass"] == f"{g1}_bias").astype(int)
    erp_df[f"flag_erpp_bias_{g2}"] = (erp_df["sexClass"] == f"{g2}_bias").astype(int)
    erp_df[f"flag_erpp_limited_{g1}"] = (erp_df["sexClass"] == f"{g1}_limited").astype(int)
    erp_df[f"flag_erpp_limited_{g2}"] = (erp_df["sexClass"] == f"{g2}_limited").astype(int)

    for status in ["w_anno_ujc", "w_anno_erp", "w_novel_erp", "w_ism"]:
        has_status = erp_df[f"flag_erpp_{status}"] == 1
        erp_df[f"flag_erpp_{status}_analyzable"] = (has_status & analyzable).astype(int)
        for group in (g1, g2):
            erp_df[f"flag_erpp_{status}_bias_{group}"] = (has_status & (erp_df["sexClass"] == f"{group}_bias")).astype(int)

    erp_df.to_csv(args.output, index=False)
    with open(args.log_file, "a") as log:
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {erp_df.shape[0]:,} rows, {erp_df.shape[1]} columns\n")

if __name__ == "__main__":
    main()