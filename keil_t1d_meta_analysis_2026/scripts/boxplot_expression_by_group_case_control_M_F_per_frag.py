#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 21:28:31 2025

@author: nkeil
"""

"""
Per-feature (exon/fragment) boxplots by group for a given gene.
Creates ONE plot per featureID instead of a single gene-wide figure.

Input CSV columns required: sampleID, featureID, <metric>
Group is parsed from sampleID suffix ...-<celltype>_<sex>_(control|case)
e.g., XXX-CD4_f_control or YYY-CD8_m_case

Examples:
  # Write one PNG per feature to a directory
  python feature_boxplots_per_fragment.py --in data.csv --gene ENSG00000000419 --metric apn --out out_plots/

  # Write a single multi-page PDF (one page per feature)
  python feature_boxplots_per_fragment.py --in data.csv --gene ENSG00000000419 --metric apn --out gene_features.pdf
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


def parse_args():
    p = argparse.ArgumentParser(
        description="One boxplot per exon/fragment for a given gene (CSV-only)."
    )
    p.add_argument("--in", dest="in_csv", required=True,
                   help="Input CSV with columns: sampleID,featureID,<metric>")
    p.add_argument("--gene", required=True,
                   help="Gene ID prefix to match (e.g., ENSG00000000419)")
    p.add_argument("--metric", default="apn",
                   help="Numeric column to plot on y-axis (default: apn)")
    p.add_argument("--out", default=None,
                   help="Output target. If endswith .pdf -> multi-page PDF. "
                        "If an existing or new directory path -> one PNG per feature. "
                        "If omitted -> show all figures interactively (careful if many).")
    p.add_argument("--min-groups", type=int, default=1,
                   help="Require at least this many groups (of the four) with data to plot a feature (default: 1).")
    return p.parse_args()


def parse_group(sample_id):
    # Expect ...-<celltype>_<sex>_(control|case) at end
    m = re.search(r"-[^_-]+_([fm])_(control|case)$", str(sample_id))
    if not m:
        return None
    sex, status = m.group(1), m.group(2)
    return f"{sex}_{status}"


def extract_gene_prefix(feature_id):
    # 'ENSG...:ER1:EF1' -> 'ENSG...'
    return str(feature_id).split(":")[0]


def feature_short_label(feature_id):
    # 'ENSG...:ER*:EF*' -> 'ER*:EF*'
    parts = str(feature_id).split(":")
    if len(parts) >= 3:
        return ":".join(parts[1:3])
    return str(feature_id)


WANTED_GROUPS = ["f_control", "f_case", "m_control", "m_case"]  # desired order
GROUP_SPECS = {
    "f_control": {"color": (0.95, 0.4, 0.4), "hatch": "",    "label": "Female control"},
    "f_case":    {"color": (0.75, 0.00, 0.00), "hatch": "//", "label": "Female case"},
    "m_control": {"color": (0.40, 0.60, 0.95), "hatch": "",    "label": "Male control"},
    "m_case":    {"color": (0.00, 0.20, 0.60), "hatch": "//", "label": "Male case"},
}


def load_and_filter(args):
    df = pd.read_csv(args.in_csv)
    required = {"sampleID", "featureID", args.metric}
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns: {missing}")

    # Gene & group parsing
    df["gene"] = df["featureID"].astype(str).map(extract_gene_prefix)
    df = df[df["gene"] == args.gene].copy()
    if df.empty:
        sys.exit(f"No rows found for gene '{args.gene}' in {args.in_csv}")

    df["group"] = df["sampleID"].map(parse_group)
    df = df[df["group"].isin(WANTED_GROUPS)].copy()
    if df.empty:
        sys.exit("No rows matched the expected groups (f_control, f_case, m_control, m_case).")

    # Ensure numeric metric
    df[args.metric] = pd.to_numeric(df[args.metric], errors="coerce")

    # Keep only features with any data
    has_any = df.groupby("featureID")[args.metric].apply(lambda s: s.notna().any())
    valid_features = has_any[has_any].index.tolist()
    if not valid_features:
        sys.exit("After filtering, no features have any data in the selected groups.")

    return df, valid_features


def plot_single_feature(ax, df_feat, metric, gene, feature_id):
    """
    df_feat: rows for one featureID, already filtered to WANTED_GROUPS
    """
    groups = []
    positions = []
    colors = []
    hatches = []
    labels = []

    x = 1
    for g in WANTED_GROUPS:
        vals = df_feat.loc[df_feat["group"] == g, metric].dropna().values
        if len(vals) == 0:
            continue
        groups.append(vals)
        positions.append(x)
        colors.append(GROUP_SPECS[g]["color"])
        hatches.append(GROUP_SPECS[g]["hatch"])
        labels.append(GROUP_SPECS[g]["label"])
        x += 1

    if len(groups) == 0:
        return False  # nothing to plot

    # Draw boxplot
    bp = ax.boxplot(
        groups,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
    )

    # Style boxes
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(colors[i])
        patch.set_edgecolor("black")
        if hatches[i]:
            patch.set_hatch(hatches[i])

    for med in bp["medians"]:
        med.set_color("black")
        med.set_linewidth(1.2)
    for whisker in bp["whiskers"]:
        whisker.set_color("black")
    for cap in bp["caps"]:
        cap.set_color("black")

    # Axis labels & ticks
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel(metric)
    ax.set_title(f"{gene} — {feature_short_label(feature_id)} ({metric})")

    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # Legend (unique)
    seen, uniq = set(), []
    for i, lbl in enumerate(labels):
        if lbl not in seen:
            seen.add(lbl)
            uniq.append(Patch(facecolor=colors[i], edgecolor="black",
                              hatch=hatches[i], label=lbl))
    ax.legend(handles=uniq, loc="best", frameon=True)

    return True


def main():
    args = parse_args()
    df, valid_features = load_and_filter(args)

    # Filter features based on how many groups have data (optional gate)
    metric = args.metric
    kept = []
    for fid in valid_features:
        grp_counts = (
            df.loc[df["featureID"] == fid]
              .groupby("group")[metric].apply(lambda s: s.notna().sum())
        )
        n_groups_with_data = int((grp_counts > 0).sum())
        if n_groups_with_data >= args.min_groups:
            kept.append(fid)

    if not kept:
        sys.exit(f"No features meet min-groups >= {args.min_groups}.")

    # Decide output mode
    if args.out is None:
        # Interactive: show a window per feature
        for fid in kept:
            fig, ax = plt.subplots(figsize=(6, 5))
            ok = plot_single_feature(ax, df[df["featureID"] == fid], metric, args.gene, fid)
            if ok:
                fig.tight_layout()
                plt.show()
            else:
                plt.close(fig)
        return

    out_path = Path(args.out)

    # Case 1: multi-page PDF
    if out_path.suffix.lower() == ".pdf":
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(out_path) as pdf:
            for fid in kept:
                fig, ax = plt.subplots(figsize=(6, 5))
                ok = plot_single_feature(ax, df[df["featureID"] == fid], metric, args.gene, fid)
                if ok:
                    fig.tight_layout()
                    pdf.savefig(fig, dpi=300, bbox_inches="tight")
                plt.close(fig)
        print(f"Wrote multi-page PDF: {out_path}")
        return

    # Case 2: directory of PNGs
    # If user passed a directory or a path without a known image suffix: treat as directory
    img_suffixes = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
    if (out_path.is_dir()) or (out_path.suffix.lower() not in img_suffixes):
        out_dir = out_path if out_path.suffix == "" or out_path.is_dir() else out_path.with_suffix("")
        out_dir.mkdir(parents=True, exist_ok=True)
        for fid in kept:
            fig, ax = plt.subplots(figsize=(6, 5))
            ok = plot_single_feature(ax, df[df["featureID"] == fid], metric, args.gene, fid)
            if ok:
                fig.tight_layout()
                short = feature_short_label(fid).replace("/", "_")
                fpath = out_dir / f"{args.gene}__{short}.png"
                fig.savefig(fpath, dpi=300, bbox_inches="tight")
            plt.close(fig)
        print(f"Wrote {len(kept)} PNG(s) to: {out_dir}")
        return

    # Case 3: single image path was given — not meaningful for multiple features.
    # We’ll write the FIRST feature to that path and put the rest next to it numbered.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    base = out_path.with_suffix("")
    ext = out_path.suffix
    for i, fid in enumerate(kept):
        fig, ax = plt.subplots(figsize=(6, 5))
        ok = plot_single_feature(ax, df[df["featureID"] == fid], metric, args.gene, fid)
        if ok:
            fig.tight_layout()
            target = base if i == 0 else Path(f"{base}_{i+1}")
            fig.savefig(str(target) + ext, dpi=300, bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote {len(kept)} image(s) starting with: {out_path}")


if __name__ == "__main__":
    main()
