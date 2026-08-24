#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate flagged UJC propReads rows to gene-level UJC summary columns."""

import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Aggregate flagged UJC propReads rows to gene-level UJC columns.")
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

    analyzable = df["flag_analyzable"] == 1
    for status in ["anno_ujc", "ujc_w_anno_erp", "ujc_w_novel_erp"]:
        df[f"flag_{status}_analyzable"] = (df[f"flag_{status}"] == 1) & analyzable
    df["flag_ism_analyzable"] = (df["flag_ism_ujc"] == 1) & analyzable

    for group in (g1, g2):
        is_bias = df["sexClass"] == f"{group}_bias"
        df[f"flag_anno_ujc_bias_{group}"] = (df["flag_anno_ujc"] == 1) & is_bias
        df[f"flag_ujc_w_anno_erp_bias_{group}"] = (df["flag_ujc_w_anno_erp"] == 1) & is_bias
        df[f"flag_ujc_w_novel_erp_bias_{group}"] = (df["flag_ujc_w_novel_erp"] == 1) & is_bias
        df[f"flag_ism_bias_{group}"] = (df["flag_ism_ujc"] == 1) & is_bias

    gene_summary = df.groupby("geneID", sort=False).agg(
        num_ujc=("jxnHash", "nunique"),
        num_anno_ujc=("flag_anno_ujc", "sum"),
        num_ujc_w_anno_erp=("flag_ujc_w_anno_erp", "sum"),
        num_ujc_w_novel_erp=("flag_ujc_w_novel_erp", "sum"),
        num_ism_ujc=("flag_ism_ujc", "sum"),
        num_MonoExon=("flag_MonoExon", "sum"),
        sumReads=("sumReads_T", "sum"),
        **{f"sumReads_{g1}": (f"sumReads_{g1}", "sum")},
        **{f"sumReads_{g2}": (f"sumReads_{g2}", "sum")},
        sumReads_MonoExon=("sumReads_monoExon", "sum"),
        num_ujc_analyzable=("flag_analyzable", "sum"),
        num_anno_ujc_analyzable=("flag_anno_ujc_analyzable", "sum"),
        num_ujc_w_anno_erp_analyzable=("flag_ujc_w_anno_erp_analyzable", "sum"),
        num_ujc_w_novel_erp_analyzable=("flag_ujc_w_novel_erp_analyzable", "sum"),
        num_ism_analyzable=("flag_ism_analyzable", "sum"),
        **{f"num_ujc_limited_{g1}": (f"flag_ujc_limited_{g1}", "sum")},
        **{f"num_ujc_limited_{g2}": (f"flag_ujc_limited_{g2}", "sum")},
        **{f"num_ujc_bias_{g1}": (f"flag_ujc_bias_{g1}", "sum")},
        **{f"num_ujc_bias_{g2}": (f"flag_ujc_bias_{g2}", "sum")},
        **{f"num_anno_ujc_bias_{g1}": (f"flag_anno_ujc_bias_{g1}", "sum")},
        **{f"num_anno_ujc_bias_{g2}": (f"flag_anno_ujc_bias_{g2}", "sum")},
        **{f"num_ujc_w_anno_erp_bias_{g1}": (f"flag_ujc_w_anno_erp_bias_{g1}", "sum")},
        **{f"num_ujc_w_anno_erp_bias_{g2}": (f"flag_ujc_w_anno_erp_bias_{g2}", "sum")},
        **{f"num_ujc_w_novel_erp_bias_{g1}": (f"flag_ujc_w_novel_erp_bias_{g1}", "sum")},
        **{f"num_ujc_w_novel_erp_bias_{g2}": (f"flag_ujc_w_novel_erp_bias_{g2}", "sum")},
        **{f"num_ism_bias_{g1}": (f"flag_ism_bias_{g1}", "sum")},
        **{f"num_ism_bias_{g2}": (f"flag_ism_bias_{g2}", "sum")},
        mstExpr_UJC=("jxnHash", lambda x: x[df.loc[x.index, "flag_mstExpr_UJC_T"] == 1].iloc[0] if (df.loc[x.index, "flag_mstExpr_UJC_T"] == 1).any() else np.nan),
        flag_mstExpr_UJC_tie=("flag_mstExpr_UJC_T", lambda x: int(x.sum() > 1)),
        **{f"mstExpr_UJC_{g1}": ("jxnHash", lambda x, g=g1: x[df.loc[x.index, f"flag_mstExpr_UJC_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_UJC_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_UJC_{g1}_tie": (f"flag_mstExpr_UJC_{g1}", lambda x: int(x.sum() > 1))},
        **{f"mstExpr_UJC_{g2}": ("jxnHash", lambda x, g=g2: x[df.loc[x.index, f"flag_mstExpr_UJC_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_UJC_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_UJC_{g2}_tie": (f"flag_mstExpr_UJC_{g2}", lambda x: int(x.sum() > 1))},
        mstExpr_annoUJC=("jxnHash", lambda x: x[df.loc[x.index, "flag_mstExpr_annoUJC_T"] == 1].iloc[0] if (df.loc[x.index, "flag_mstExpr_annoUJC_T"] == 1).any() else np.nan),
        flag_mstExpr_annoUJC_T_tie=("flag_mstExpr_annoUJC_T", lambda x: int(x.sum() > 1)),
        **{f"mstExpr_annoUJC_{g1}": ("jxnHash", lambda x, g=g1: x[df.loc[x.index, f"flag_mstExpr_annoUJC_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_annoUJC_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_annoUJC_{g1}_tie": (f"flag_mstExpr_annoUJC_{g1}", lambda x: int(x.sum() > 1))},
        **{f"mstExpr_annoUJC_{g2}": ("jxnHash", lambda x, g=g2: x[df.loc[x.index, f"flag_mstExpr_annoUJC_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_annoUJC_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_annoUJC_{g2}_tie": (f"flag_mstExpr_annoUJC_{g2}", lambda x: int(x.sum() > 1))},
        propReads_mstExpr_UJC=("T_read_proportion", "max"),
        **{f"propReads_mstExpr_UJC_{g1}": (f"{g1}_read_proportion", "max")},
        **{f"propReads_mstExpr_UJC_{g2}": (f"{g2}_read_proportion", "max")},
        propReads_mstExpr_annoUJC=("T_read_proportion_anno", "max"),
        **{f"propReads_mstExpr_annoUJC_{g1}": (f"{g1}_read_proportion_anno", "max")},
        **{f"propReads_mstExpr_annoUJC_{g2}": (f"{g2}_read_proportion_anno", "max")},
    ).reset_index()

    gene_summary.to_csv(args.output, index=False)

    with open(args.log_file, "w") as log:
        log.write(f"Aggregate UJC propReads to gene summary: {args.species}\n")
        log.write(f"Input: {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {gene_summary.shape[0]:,} genes, {gene_summary.shape[1]} columns\n")

if __name__ == "__main__":
    main()