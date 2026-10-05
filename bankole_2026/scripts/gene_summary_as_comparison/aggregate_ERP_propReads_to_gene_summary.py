#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate flagged ERP_plus propReads rows to gene-level ERP summary columns."""

import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Aggregate flagged ERP propReads rows to gene-level ERP columns.")
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

    agg = dict(
        num_ERP=("ERP", "nunique"),
        num_ERPp=("ERP_plus", "nunique"),
        num_ERPp_w_anno_ujc=("flag_erpp_w_anno_ujc", "sum"),
        num_ERPp_w_anno_erp=("flag_erpp_w_anno_erp", "sum"),
        num_ERPp_w_novel_erp=("flag_erpp_w_novel_erp", "sum"),
        num_ERPp_w_ism=("flag_erpp_w_ism", "sum"),
        numExon_GM=("numExon_GM", "first"),
        mstExpr_ERPp=("ERP_plus", lambda x: x[df.loc[x.index, "flag_mstExpr_ERPp_T"] == 1].iloc[0] if (df.loc[x.index, "flag_mstExpr_ERPp_T"] == 1).any() else np.nan),
        flag_mstExpr_ERPp_tie=("flag_mstExpr_ERPp_T", lambda x: int(x.sum() > 1)),
        **{f"mstExpr_ERPp_{g1}": ("ERP_plus", lambda x, g=g1: x[df.loc[x.index, f"flag_mstExpr_ERPp_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_ERPp_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_ERPp_{g1}_tie": (f"flag_mstExpr_ERPp_{g1}", lambda x: int(x.sum() > 1))},
        **{f"mstExpr_ERPp_{g2}": ("ERP_plus", lambda x, g=g2: x[df.loc[x.index, f"flag_mstExpr_ERPp_{g}"] == 1].iloc[0] if (df.loc[x.index, f"flag_mstExpr_ERPp_{g}"] == 1).any() else np.nan)},
        **{f"flag_mstExpr_ERPp_{g2}_tie": (f"flag_mstExpr_ERPp_{g2}", lambda x: int(x.sum() > 1))},
        propReads_mstExpr_ERPp=("T_read_proportion", "max"),
        **{f"propReads_mstExpr_ERPp_{g1}": (f"{g1}_read_proportion", "max")},
        **{f"propReads_mstExpr_ERPp_{g2}": (f"{g2}_read_proportion", "max")},
        num_ERPp_analyzable=("flag_analyzable", "sum"),
        **{f"num_ERPp_bias_{g1}": (f"flag_erpp_bias_{g1}", "sum")},
        **{f"num_ERPp_bias_{g2}": (f"flag_erpp_bias_{g2}", "sum")},
        num_ERPp_w_anno_ujc_analyzable=("flag_erpp_w_anno_ujc_analyzable", "sum"),
        num_ERPp_w_anno_erp_analyzable=("flag_erpp_w_anno_erp_analyzable", "sum"),
        num_ERPp_w_novel_erp_analyzable=("flag_erpp_w_novel_erp_analyzable", "sum"),
        num_ERPp_w_ism_analyzable=("flag_erpp_w_ism_analyzable", "sum"),
        **{f"num_ERPp_w_anno_ujc_bias_{g1}": (f"flag_erpp_w_anno_ujc_bias_{g1}", "sum")},
        **{f"num_ERPp_w_anno_ujc_bias_{g2}": (f"flag_erpp_w_anno_ujc_bias_{g2}", "sum")},
        **{f"num_ERPp_w_anno_erp_bias_{g1}": (f"flag_erpp_w_anno_erp_bias_{g1}", "sum")},
        **{f"num_ERPp_w_anno_erp_bias_{g2}": (f"flag_erpp_w_anno_erp_bias_{g2}", "sum")},
        **{f"num_ERPp_w_novel_erp_bias_{g1}": (f"flag_erpp_w_novel_erp_bias_{g1}", "sum")},
        **{f"num_ERPp_w_novel_erp_bias_{g2}": (f"flag_erpp_w_novel_erp_bias_{g2}", "sum")},
        **{f"num_ERPp_w_ism_bias_{g1}": (f"flag_erpp_w_ism_bias_{g1}", "sum")},
        **{f"num_ERPp_w_ism_bias_{g2}": (f"flag_erpp_w_ism_bias_{g2}", "sum")},
        **{f"num_ERPp_limited_{g1}": (f"flag_erpp_limited_{g1}", "sum")},
        **{f"num_ERPp_limited_{g2}": (f"flag_erpp_limited_{g2}", "sum")},
    )

    gene_summary = df.groupby("geneID", sort=False).agg(**agg).reset_index()
    gene_summary.to_csv(args.output, index=False)

    with open(args.log_file, "w") as log:
        log.write(f"Aggregate ERP propReads to gene summary: {args.species}\n")
        log.write(f"Input: {args.input}\n")
        log.write(f"Output: {args.output}\n")
        log.write(f"Final: {gene_summary.shape[0]:,} genes, {gene_summary.shape[1]} columns\n")

if __name__ == "__main__":
    main()