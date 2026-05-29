#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-fragment (featureID) boxplots for two groups: cases (black) vs controls (grey).

- Input CSV columns: sampleID, featureID, <metric>
- Filters to a single gene (by featureID prefix before the first ':')
- No positions/ordering file is supported here (fragments are ordered alphanumerically)
- Keeps only fragments with data in at least one of the two groups
- Each OUTPUT FIGURE contains ONLY ONE fragment (featureID)
- X-axis shows two boxes (case, control) for that fragment (only groups with data are drawn)

Examples:
  # Write one PNG per fragment to a directory
  python per_fragment_cases_controls.py \
    --in data.csv --gene ENSG00000000419 --metric apn --out out_plots/

  # Write a single multi-page PDF (one page per fragment)
  python per_fragment_cases_controls.py \
    --in data.csv --gene ENSG00000000419 --metric apn --out gene_fragments.pdf
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Use a non-GUI backend to avoid XInput warnings and headless errors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch


def parse_args():
    p = argparse.ArgumentParser(
        description="One boxplot per fragment (featureID) for cases vs controls (no positions file)."
    )
    p.add_argument("--in", dest="in_csv", required=True,
                   help="Input CSV with columns: sampleID,featureID,<metric>")
    p.add_argument("--gene", required=True,
                   help="Gene ID prefix to match (e.g., ENSG00000000419)")
    p.add_argument("--metric", default="apn",
                   help="Numeric column to plot on y-axis (default: apn)")
    p.add_argument("--out", default=None,
                   help="If endswith .pdf -> multi-page PDF. "
                        "If a directory -> one PNG per fragment. "
                        "If omitted -> shows interactive windows (not recommended on headless systems).")
    p.add_argument("--min-groups", type=int, default=1,
                   help="Require at least this many groups (of the two) with data to plot a fragment (default: 1).")
    p.add_argument("--dpi", type=int, default=300,
                   help="Output image DPI (default: 300).")
    return p.parse_args()


def parse_group(sample_id):
    """
    Map sampleID -> 'case' or 'control' by looking for a trailing _case/_control.
    Matches e.g.: ...-CD8_control, ...-CD8_f_case
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


def feature_short_label(feature_id):
    # 'ENSG...:ER*:EF*' -> 'ER*:EF*'
    parts = str(feature_id).split(":")
    if len(parts) >= 3:
        return ":".join(parts[1:3])
    return str(feature_id)


WANTED_GROUPS = ["case", "control"]  # order: cases (black) then controls (grey)
GROUP_SPECS = {
    "case":    {"color": "black", "label": "Case"},
    "control": {"color": "grey",  "label": "Control"},
}


def load_and_filter(args):
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

    # Parse 2-group label
    df["group"] = df["sampleID"].map(parse_group)
    df = df[df["group"].isin(WANTED_GROUPS)].copy()
    if df.empty:
        sys.exit("No rows matched the expected groups (case, control). Check sampleID naming.")

    # Ensure metric numeric
    df[args.metric] = pd.to_numeric(df[args.metric], errors="coerce")

    # Keep only fragments with any data in either group
    has_any = df.groupby("featureID")[args.metric].apply(lambda s: s.notna().any())
    valid_features = has_any[has_any].index.tolist()
    if not valid_features:
        sys.exit("After filtering, no fragments have any data in the selected groups.")

    return df, valid_features


def plot_single_fragment(ax, df_frag, metric, gene, feature_id):
    """
    df_frag: rows for one featureID, already filtered to desired groups
    """
    groups = []
    positions = []
    colors = []
    labels = []

    x = 1
    for g in WANTED_GROUPS:
        vals = df_frag.loc[df_frag["group"] == g, metric].dropna().values
        if len(vals) == 0:
            continue
        groups.append(vals)
        positions.append(x)
        colors.append(GROUP_SPECS[g]["color"])
        labels.append(GROUP_SPECS[g]["label"])
        x += 1

    if len(groups) == 0:
        return False  # nothing to plot

    bp = ax.boxplot(
        groups,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        manage_ticks=False,
        showfliers=False,
    )

    # Style
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(colors[i])
        patch.set_edgecolor("black")

    for i, med in enumerate(bp["medians"]):
        med.set_color("white" if colors[i] == "black" else "black")
        med.set_linewidth(1.2)

    for whisker in bp["whiskers"]:
        whisker.set_color("black")
    for cap in bp["caps"]:
        cap.set_color("black")

    # Axis labels & ticks
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel(metric)
    ax.set_title(f"{gene} — {feature_short_label(feature_id)} ({metric})")
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # Legend (unique)
    seen, uniq = set(), []
    for i, lbl in enumerate(labels):
        if lbl not in seen:
            seen.add(lbl)
            uniq.append(Patch(facecolor=colors[i], edgecolor="black", label=lbl))
    ax.legend(handles=uniq, loc="best", frameon=True)

    return True


def main():
    args = parse_args()
    df, valid_features = load_and_filter(args)

    # Require at least k of the two groups to have >=1 datum for the fragment
    kept = []
    metric = args.metric
    for fid in valid_features:
        grp_counts = (
            df.loc[df["featureID"] == fid]
              .groupby("group")[metric].apply(lambda s: s.notna().sum())
        )
        n_groups_with_data = int((grp_counts > 0).sum())
        if n_groups_with_data >= args.min_groups:
            kept.append(fid)

    if not kept:
        sys.exit(f"No fragments meet min-groups >= {args.min_groups}.")

    # Fragment page order: simple alphanumeric by featureID (no positions supported)
    feature_order = sorted(kept, key=str)

    # Output handling
    if args.out is None:
        # Interactive windows (be careful on headless systems)
        for fid in feature_order:
            fig, ax = plt.subplots(figsize=(6, 5))
            ok = plot_single_fragment(ax, df[df["featureID"] == fid], metric, args.gene, fid)
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
            for fid in feature_order:
                fig, ax = plt.subplots(figsize=(6, 5))
                ok = plot_single_fragment(ax, df[df["featureID"] == fid], metric, args.gene, fid)
                if ok:
                    fig.tight_layout()
                    pdf.savefig(fig, dpi=args.dpi, bbox_inches="tight")
                plt.close(fig)
        print(f"Wrote multi-page PDF: {out_path}")
        return

    # Case 2: directory of PNGs (if path is a dir OR not a known single-image suffix)
    img_suffixes = {".png", ".jpg", ".jpeg", ".tif", ".tiff"}
    if (out_path.is_dir()) or (out_path.suffix.lower() not in img_suffixes):
        out_dir = out_path if out_path.suffix == "" or out_path.is_dir() else out_path.with_suffix("")
        out_dir.mkdir(parents=True, exist_ok=True)
        for fid in feature_order:
            fig, ax = plt.subplots(figsize=(6, 5))
            ok = plot_single_fragment(ax, df[df["featureID"] == fid], metric, args.gene, fid)
            if ok:
                fig.tight_layout()
                short = feature_short_label(fid).replace("/", "_").replace(":", "_")
                fpath = out_dir / f"{args.gene}__{short}.png"
                fig.savefig(fpath, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
        print(f"Wrote {len(feature_order)} PNG(s) to: {out_dir}")
        return

    # Case 3: single image path given — write first fragment to that path, rest numbered
    out_path.parent.mkdir(parents=True, exist_ok=True)
    base = out_path.with_suffix("")
    ext = out_path.suffix
    for i, fid in enumerate(feature_order):
        fig, ax = plt.subplots(figsize=(6, 5))
        ok = plot_single_fragment(ax, df[df["featureID"] == fid], metric, args.gene, fid)
        if ok:
            fig.tight_layout()
            target = base if i == 0 else Path(f"{base}_{i+1}")
            fig.savefig(str(target) + ext, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
    print(f"Wrote {len(feature_order)} image(s) starting with: {out_path}")


if __name__ == "__main__":
    main()
