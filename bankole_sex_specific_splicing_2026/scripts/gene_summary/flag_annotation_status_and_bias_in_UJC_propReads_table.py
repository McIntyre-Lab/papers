#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Add annotation-partition, ISM, mono-exon, and sex flags to a UJC propReads table.

Annotation counts follow the column-definition partition on the flags
(flag_anno / flag_ERPanno); ISM is a separate overlapping flag. annotation_status
is the hierarchical string (anno_UJC > ISM > anno_ERP > novel_ERP) used only by
the downstream AS-event detection to exclude ISM UJC.
"""

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Flag annotation status and bias in a UJC propReads table.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("-l", "--log-file", required=True)
    parser.add_argument("--ujc-annotation-flag", default="flag_jxnHash_in_fiveSpecies_full_anno")
    parser.add_argument("--erp-annotation-flag", default="flag_ERPanno")
    parser.add_argument("--group-class-col", default="sexClass")
    args = parser.parse_args()

    group1 = args.group1
    group2 = args.group2
    anno_flag = args.ujc_annotation_flag
    erp_flag = args.erp_annotation_flag

    df = pd.read_csv(args.input, low_memory=False)
    df = df.rename(columns={args.group_class_col: "sexClass"})

    is_anno = df[anno_flag] == 1
    is_erp_anno = df[erp_flag] == 1
    is_ism = df["flag_is_incomplete-splice_match"] == 1

    # Mutually exclusive annotation status assigned hierarchically: anno_UJC > ISM > anno_ERP > novel_ERP),
    df["annotation_status"] = np.select(
        [is_anno, is_ism, is_erp_anno],
        ["anno_UJC", "ISM", "anno_ERP"],
        default="novel_ERP",
    )

    # Count flags derived from the exclusive status: every UJC lands in exactly one bucket.
    df["flag_anno_ujc"] = (df["annotation_status"] == "anno_UJC").astype(int)
    df["flag_ujc_w_anno_erp"] = (df["annotation_status"] == "anno_ERP").astype(int)
    df["flag_ujc_w_novel_erp"] = (df["annotation_status"] == "novel_ERP").astype(int)
    df["flag_ism_ujc"] = (df["annotation_status"] == "ISM").astype(int)
    df["flag_MonoExon"] = (df["numExon_ERP"] == 1).astype(int)

    df[f"flag_ujc_bias_{group1}"] = (df["sexClass"] == f"{group1}_bias").astype(int)
    df[f"flag_ujc_bias_{group2}"] = (df["sexClass"] == f"{group2}_bias").astype(int)
    df[f"flag_ujc_limited_{group1}"] = (df["sexClass"] == f"{group1}_limited").astype(int)
    df[f"flag_ujc_limited_{group2}"] = (df["sexClass"] == f"{group2}_limited").astype(int)

    df.to_csv(args.output, index=False)

    with open(args.log_file, "w") as log:
        log.write(f"Flag UJC propReads table: {args.species}\n")
        log.write(f"Input: {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {df.shape[0]:,} rows, {df.shape[1]} columns\n")


if __name__ == "__main__":
    main()
