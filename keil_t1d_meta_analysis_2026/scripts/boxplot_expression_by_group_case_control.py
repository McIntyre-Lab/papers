#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 13:55:32 2025

@author: nkeil
"""

"""
Per-feature boxplots for two groups: cases (black) vs controls (grey).

- Input CSV columns: sampleID, featureID, <metric>
- Filters to a single gene (by featureID prefix before the first ':')
- Optionally sorts features by 'start' from a positions CSV
- Keeps only features with data in at least one of the two groups
- X-axis labels show ER*:EF* (gene ID removed)
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def parse_args():
    p = argparse.ArgumentParser(description="Per-feature boxplots for cases vs controls (CSV-only).")
    p.add_argument("--in", dest="in_csv", required=True,
                   help="Input CSV with columns: sampleID,featureID,<metric>")
    p.add_argument("--gene", required=True,
                   help="Gene ID prefix to match (e.g., ENSG00000000419)")
    p.add_argument("--metric", default="apn",
                   help="Numeric column to plot on y-axis (default: apn)")
    p.add_argument("--positions", default=None,
                   help="Optional CSV with columns: featureID,start (used to sort x-axis)")
    p.add_argument("--out", default=None,
                   help="Output image (e.g., .png/.pdf). If omitted, shows window.")
    return p.parse_args()

def parse_group(sample_id):
    """
    Map sampleID -> 'case' or 'control' by looking for a trailing _case/_control.
    Examples that match:
      ...-CD8_control
      ...-CD8_f_case
      ...-whatever_anything_control
    """
    s = str(sample_id)
    m = re.search(r"_(control|case)$", s)
    if m:
        return m.group(1)
    # Fallback: also accept '-..._(control|case)' at the very end (rare naming)
    m2 = re.search(r"(?:^|[-_])(control|case)$", s)
    return m2.group(1) if m2 else None

def extract_gene_prefix(feature_id):
    # 'ENSG...:ER1:EF1' -> 'ENSG...'
    return str(feature_id).split(":")[0]

def feature_label(feature_id):
    # Show only ER*:EF* on x-axis
    parts = str(feature_id).split(":")
    if len(parts) >= 3:
        return ":".join(parts[1:3])
    return str(feature_id)

def main():
    args = parse_args()
    df = pd.read_csv(args.in_csv)

    required = {"sampleID", "featureID", args.metric}
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"Missing required columns: {missing}")

    # Filter to selected gene
    df["gene"] = df["featureID"].astype(str).map(extract_gene_prefix)
    df = df[df["gene"] == args.gene].copy()
    if df.empty:
        sys.exit(f"No rows found for gene '{args.gene}' in {args.in_csv}")

    # Parse 2-group label: 'case' or 'control'
    df["group"] = df["sampleID"].map(parse_group)
    wanted_groups = ["case", "control"]  # desired order: cases first (black), controls second (grey)
    df = df[df["group"].isin(wanted_groups)].copy()
    if df.empty:
        sys.exit("No rows matched the expected groups (case, control). Check sampleID naming.")

    # Ensure metric numeric
    df[args.metric] = pd.to_numeric(df[args.metric], errors="coerce")

    # Keep only features that have data in at least one of the two groups
    has_any_data = df.groupby("featureID")[args.metric].apply(lambda s: s.notna().any())
    valid_features = set(has_any_data[has_any_data].index)
    if not valid_features:
        sys.exit("After filtering, no features have any data in the selected groups.")

    df = df[df["featureID"].isin(valid_features)].copy()

    # Feature order (optionally by start)
    if args.positions:
        pos = pd.read_csv(args.positions)
        if not {"featureID", "start"}.issubset(pos.columns):
            print("Positions CSV must contain 'featureID' and 'start'. Using alphanumeric sort.", file=sys.stderr)
            feature_order = sorted(valid_features, key=str)
        else:
            pos2 = pos[["featureID", "start"]].copy()
            pos2["start"] = pd.to_numeric(pos2["start"], errors="coerce")
            df = df.merge(pos2, on="featureID", how="left")
            starts = (
                df.groupby("featureID")["start"].min()
                  .reindex(sorted(valid_features, key=str))
            )
            with_start = starts[starts.notna()].sort_values(kind="mergesort").index.tolist()
            without_start = starts[starts.isna()].index.tolist()
            feature_order = with_start + without_start
    else:
        feature_order = sorted(valid_features, key=str)

    # Visual encoding: cases = black, controls = grey
    group_specs = {
        "case":    {"color": "black", "label": "Case"},
        "control": {"color": "grey",  "label": "Control"},
    }

    n_features = len(feature_order)
    base_x = np.arange(n_features, dtype=float) + 1.0
    offsets = {"case": -0.12, "control": +0.12}
    width = 0.20

    fig, ax = plt.subplots(figsize=(max(8, n_features * 0.6), 6))
    legend_handles = []

    for g in wanted_groups:
        xs, groups = [], []
        for i, feat in enumerate(feature_order):
            vals = df.loc[(df["featureID"] == feat) & (df["group"] == g), args.metric].dropna().values
            if len(vals) == 0:
                vals = np.array([np.nan])  # keep slot alignment
            xs.append(base_x[i] + offsets[g])
            groups.append(vals)

        bp = ax.boxplot(
            groups,
            positions=xs,
            widths=width,
            patch_artist=True,
            manage_ticks=False,
            showfliers=False,
        )

        col = group_specs[g]["color"]
        for patch in bp["boxes"]:
            patch.set_facecolor(col)
            patch.set_edgecolor("black")
        for median in bp["medians"]:
            median.set_color("white" if col == "black" else "black")
            median.set_linewidth(1.2)
        for whisker in bp["whiskers"]:
            whisker.set_color("black")
        for cap in bp["caps"]:
            cap.set_color("black")

        legend_handles.append(Patch(facecolor=col, edgecolor="black", label=group_specs[g]["label"]))

    ax.set_xticks(base_x)
    ax.set_xticklabels([feature_label(f) for f in feature_order], rotation=45, ha="right")
    ax.set_xlabel("Feature ID" + (" (sorted by start)" if args.positions else ""))
    ax.set_ylabel(args.metric)
    ax.set_title(f"{args.gene} — per-feature boxplots (cases vs controls, {args.metric})")

    # Legend without a title and with unique labels
    seen, uniq = set(), []
    for h in legend_handles:
        if h.get_label() not in seen:
            seen.add(h.get_label())
            uniq.append(h)
    ax.legend(handles=uniq, loc="best", frameon=True)

    ax.grid(axis="y", linestyle=":", alpha=0.4)
    fig.tight_layout()

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Wrote {out_path}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
