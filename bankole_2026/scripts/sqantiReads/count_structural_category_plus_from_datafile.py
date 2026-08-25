#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 15:33:40 2025

@author: nkeil
"""
# build_per_sample_struct_tables.py
import os
import pandas as pd

# Paths
ind  = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis"
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/sqanti_reads_analysis"

# Species tags present
species_tags = {
    "mel": "dmel6",
    "sim": "dsim2",
    "yak": "dyak2",
    "san": "dsan1",
    "ser": "dser1",
}

cnts_prefix = "rawCnts"

# Preferred order for output columns
preferred_order = [
    "FSM", "ISM", "NIC", "NNC_annoERP", "NNC_no_annoERP",
    "AS", "FUS", "GENIC", "GI", "INTER", "AMBIG"
]

#Function to count structural category plus for each sample
def sample_category_read_sums(df, prefix=cnts_prefix):
    """Per-sample SUMS of raw counts grouped by structural_category_plus."""
    if "structural_category_plus" not in df.columns:
        raise ValueError("Missing 'structural_category_plus' in input CSV.")
    cnt_cols = [c for c in df.columns if c.startswith(prefix)]
    if not cnt_cols:
        raise ValueError(f"No columns starting with '{prefix}' found.")

    # numeric counts
    df[cnt_cols] = df[cnt_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    # column -> sample name (strip prefix)
    col_to_sample = {c: c.split(prefix, 1)[1].lstrip("_") for c in cnt_cols}

    # sum reads per category for each sample
    cat = df["structural_category_plus"]
    rows = {}
    for c in cnt_cols:
        sums = df[c].groupby(cat).sum()
        rows[col_to_sample[c]] = sums.to_dict()

    out = pd.DataFrame.from_dict(rows, orient="index").fillna(0)

    # order columns, add total
    existing = list(out.columns)
    ordered = [x for x in preferred_order if x in existing] + [x for x in existing if x not in preferred_order]
    out = out[ordered]
    out["total"] = out.sum(axis=1)
    return out

#Function to add percent structural category columns
def add_percent_columns(df_with_counts):
    """Add perc_<CAT> columns to a counts table that already has 'total'."""
    if "total" not in df_with_counts.columns:
        raise ValueError("'total' column missing.")
    counts_cols = [c for c in df_with_counts.columns if c != "total" and not c.startswith("perc_")]
    for cat in counts_cols:
        df_with_counts[f"perc_{cat}"] = (df_with_counts[cat] / df_with_counts["total"]).where(df_with_counts["total"] > 0, 0) * 100.0
    return df_with_counts

# Make per species and combined tables
tables = []
for _, tag in species_tags.items():
    in_file = os.path.join(ind, f"datafile_jxnHash_{tag}_w_struct_category_plus.csv")
    if not os.path.exists(in_file):
        continue
    df = pd.read_csv(in_file, low_memory=False)
    per_sp = sample_category_read_sums(df)
    per_sp = add_percent_columns(per_sp)
    #per_sp.to_csv(os.path.join(outd, f"per_sample_struct_readSums_percent_{tag}.csv"))
    tables.append(per_sp)

if not tables:
    raise SystemExit("No species tables were built—check your input files.")

combined = pd.concat(tables, axis=0).fillna(0)

# Re-enforce counts-column order & keep 'total' at end of counts block.
count_cols_existing = [c for c in combined.columns if not c.startswith("perc_") and c != "total"]
ordered_counts = [x for x in preferred_order if x in count_cols_existing] + [x for x in count_cols_existing if x not in preferred_order]
combined = combined[ordered_counts + ["total"] + [c for c in combined.columns if c.startswith("perc_")]]

combined_out = os.path.join(outd, "structural_category_plus_counts_by_sample.csv")
combined.to_csv(combined_out)

print(f"✓ Wrote (per-species): per_sample_struct_readSums_percent_{tag}.csv")
print(f"✓ Wrote (combined): {combined_out}")