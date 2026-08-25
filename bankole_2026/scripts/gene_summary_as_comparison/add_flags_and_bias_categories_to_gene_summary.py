#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Add flags and sex-bias category columns to a gene summary table.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("-l", "--log-file", required=True)
    args = parser.parse_args()

    g1 = args.group1
    g2 = args.group2

    df = pd.read_csv(args.input, low_memory=False)

    df["flag_gff5_gene"] = df["geneID"].astype(str).str.startswith("GFF5").astype(int)
    df["flag_only_ISM"] = (df["num_ERPp"] == df["num_ERPp_w_ism"]).astype(int)

    # sexBias: any biased UJC or ERP_plus in each direction (>= 1)
    f_bias = (df[f"num_ujc_bias_{g1}"] >= 1) | (df[f"num_ERPp_bias_{g1}"] >= 1)
    m_bias = (df[f"num_ujc_bias_{g2}"] >= 1) | (df[f"num_ERPp_bias_{g2}"] >= 1)
    analyzable = (df["num_ujc_analyzable"] >= 1) | (df["num_ERPp_analyzable"] >= 1)
    df["sexBias"] = np.select(
        [f_bias & m_bias, f_bias & ~m_bias, m_bias & ~f_bias, analyzable],
        ["B", g1, g2, "unbiased"], default="not_evaluated",
    )

    # sexBias_annoUJC: direction from annotated UJC only
    fa = df[f"num_anno_ujc_bias_{g1}"] >= 1
    ma = df[f"num_anno_ujc_bias_{g2}"] >= 1
    ana = df["num_anno_ujc_analyzable"] >= 1
    df["sexBias_annoUJC"] = np.select(
        [fa & ma, fa & ~ma, ma & ~fa, ana],
        ["B", g1, g2, "unbiased"], default="not_evaluated",
    )

    # sexBias_anno_status: most-annotated bias tier present (tiers match the existing pipeline;
    # note num_ERPp_w_anno_ujc bias is grouped under the anno_ERP tier, as in the published version).
    anno_ujc_bias = (df[f"num_anno_ujc_bias_{g1}"] >= 1) | (df[f"num_anno_ujc_bias_{g2}"] >= 1)
    anno_erp_bias = (df[f"num_ERPp_w_anno_ujc_bias_{g1}"] >= 1) | (df[f"num_ERPp_w_anno_ujc_bias_{g2}"] >= 1) | \
                    (df[f"num_ERPp_w_anno_erp_bias_{g1}"] >= 1) | (df[f"num_ERPp_w_anno_erp_bias_{g2}"] >= 1) | \
                    (df[f"num_ujc_w_anno_erp_bias_{g1}"] >= 1) | (df[f"num_ujc_w_anno_erp_bias_{g2}"] >= 1)
    novel_erp_bias = (df[f"num_ERPp_w_novel_erp_bias_{g1}"] >= 1) | (df[f"num_ERPp_w_novel_erp_bias_{g2}"] >= 1) | \
                     (df[f"num_ujc_w_novel_erp_bias_{g1}"] >= 1) | (df[f"num_ujc_w_novel_erp_bias_{g2}"] >= 1)
    ism_bias = (df[f"num_ism_bias_{g1}"] >= 1) | (df[f"num_ism_bias_{g2}"] >= 1) | \
               (df[f"num_ERPp_w_ism_bias_{g1}"] >= 1) | (df[f"num_ERPp_w_ism_bias_{g2}"] >= 1)
    any_analyzable = (df["num_anno_ujc_analyzable"] >= 1) | (df["num_ujc_w_anno_erp_analyzable"] >= 1) | \
                     (df["num_ujc_w_novel_erp_analyzable"] >= 1) | (df["num_ism_analyzable"] >= 1) | \
                     (df["num_ERPp_w_anno_ujc_analyzable"] >= 1) | (df["num_ERPp_w_anno_erp_analyzable"] >= 1) | \
                     (df["num_ERPp_w_novel_erp_analyzable"] >= 1) | (df["num_ERPp_w_ism_analyzable"] >= 1)
    df["sexBias_anno_status"] = np.select(
        [anno_ujc_bias, anno_erp_bias, novel_erp_bias, ism_bias, any_analyzable],
        ["anno_UJC_bias", "anno_ERP_bias", "novel_ERP_bias", "ISM_bias", "unbiased"],
        default="not_evaluated",
    )

    ordered_columns = [
        "geneID",
        "num_ujc", "num_anno_ujc", "num_ujc_w_anno_erp", "num_ujc_w_novel_erp", "num_ism_ujc", "num_MonoExon",
        "sumReads", f"sumReads_{g1}", f"sumReads_{g2}", "sumReads_MonoExon",
        "num_ujc_analyzable", "num_anno_ujc_analyzable", "num_ujc_w_anno_erp_analyzable",
        "num_ujc_w_novel_erp_analyzable", "num_ism_analyzable",
        f"num_ujc_limited_{g1}", f"num_ujc_limited_{g2}", f"num_ujc_bias_{g1}", f"num_ujc_bias_{g2}",
        f"num_anno_ujc_bias_{g1}", f"num_anno_ujc_bias_{g2}",
        f"num_ujc_w_anno_erp_bias_{g1}", f"num_ujc_w_anno_erp_bias_{g2}",
        f"num_ujc_w_novel_erp_bias_{g1}", f"num_ujc_w_novel_erp_bias_{g2}",
        f"num_ism_bias_{g1}", f"num_ism_bias_{g2}",
        "mstExpr_UJC", f"mstExpr_UJC_{g1}", f"mstExpr_UJC_{g2}",
        "mstExpr_annoUJC", f"mstExpr_annoUJC_{g1}", f"mstExpr_annoUJC_{g2}",
        "ERPp_for_mstExpr_UJC", f"ERPp_for_mstExpr_UJC_{g1}", f"ERPp_for_mstExpr_UJC_{g2}",
        "ERPp_for_mstExpr_annoUJC", f"ERPp_for_mstExpr_annoUJC_{g1}", f"ERPp_for_mstExpr_annoUJC_{g2}",
        "propReads_mstExpr_UJC", f"propReads_mstExpr_UJC_{g1}", f"propReads_mstExpr_UJC_{g2}",
        "propReads_mstExpr_annoUJC", f"propReads_mstExpr_annoUJC_{g1}", f"propReads_mstExpr_annoUJC_{g2}",
        "num_ERP", "num_ERPp", "num_ERPp_w_anno_ujc", "num_ERPp_w_anno_erp", "num_ERPp_w_novel_erp",
        "num_ERPp_w_ism", "numExon_GM",
        "mstExpr_ERPp", f"mstExpr_ERPp_{g1}", f"mstExpr_ERPp_{g2}",
        "propReads_mstExpr_ERPp", f"propReads_mstExpr_ERPp_{g1}", f"propReads_mstExpr_ERPp_{g2}",
        "num_ERPp_analyzable", f"num_ERPp_bias_{g1}", f"num_ERPp_bias_{g2}",
        "num_ERPp_w_anno_ujc_analyzable", "num_ERPp_w_anno_erp_analyzable",
        "num_ERPp_w_novel_erp_analyzable", "num_ERPp_w_ism_analyzable",
        f"num_ERPp_w_anno_ujc_bias_{g1}", f"num_ERPp_w_anno_ujc_bias_{g2}",
        f"num_ERPp_w_anno_erp_bias_{g1}", f"num_ERPp_w_anno_erp_bias_{g2}",
        f"num_ERPp_w_novel_erp_bias_{g1}", f"num_ERPp_w_novel_erp_bias_{g2}",
        f"num_ERPp_w_ism_bias_{g1}", f"num_ERPp_w_ism_bias_{g2}",
        f"num_ERPp_limited_{g1}", f"num_ERPp_limited_{g2}",
        "flag_gff5_gene", "flag_only_ISM", "sexBias", "sexBias_anno_status", "sexBias_annoUJC",
        "flag_mstExpr_ERPp_tie", f"flag_mstExpr_ERPp_{g1}_tie", f"flag_mstExpr_ERPp_{g1}_tie",
        "flag_mstExpr_UJC_tie", f"flag_mstExpr_UJC_{g1}_tie", f"flag_mstExpr_UJC_{g2}_tie",
        "flag_mstExpr_annoUJC_tie", f"flag_mstExpr_annoUJC_{g1}_tie", f"flag_mstExpr_annoUJC_{g2}_tie"
    ]

    df = df[ordered_columns]
    df.to_csv(args.output, index=False)

    with open(args.log_file, "w") as log:
        log.write(f"Add flags and bias categories: {args.species}\n")
        log.write(f"Input: {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {df.shape[0]:,} genes, {df.shape[1]} columns\n")
        log.write(f"sexBias:\n{df['sexBias'].value_counts(dropna=False)}\n")
        log.write(f"sexBias_anno_status:\n{df['sexBias_anno_status'].value_counts(dropna=False)}\n")


if __name__ == "__main__":
    main()
