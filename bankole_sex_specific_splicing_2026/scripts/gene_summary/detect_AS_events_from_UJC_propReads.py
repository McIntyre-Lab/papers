#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""All-to-all AS event detection by gene, collapsed to gene-level AS columns.

For every unique pair of non-ISM UJC in a gene, pairwise_as_detection() flags
structural differences (alt 5p / 3p / ERSkip / donor-acceptor / IR / dataOnly).
Each event's gene-level value is the best (most annotated) tier observed across
event-positive pairs: anno_UJC > anno_ERP > novel_ERP, else 'none'. A parallel
set of *_FvM columns restricts to group1-biased vs group2-biased pairs.
Output is a standalone AS table (not merged into the gene summary).
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import pandas as pd
from pairwise_as_detection import pairwise_as_detection, erp_to_set


# pairwise_as_detection() returns these 0/1 flag keys; each maps to an output
# column name in the AS table.
EVENT_OUTPUT = {
    "flag_alt_dataOnlyExon": "alt_dataOnly",
    "flag_alt_IR": "alt_IR",
    "flag_alt_donorAcceptor": "alt_donor_acceptor",
    "flag_alt_5pER": "alt_5p_ER",
    "flag_alt_3pER": "alt_3p_ER",
    "flag_alt_skipER": "alt_ERSkip",

}

EVENT_TYPES = list(EVENT_OUTPUT.keys())

ANNOTATION_TIER_MAP = {"anno_UJC": 1, "anno_ERP": 2, "novel_ERP": 3}
TIER_OUTPUT = {1: "anno_UJC", 2: "anno_ERP", 3: "novel_ERP"}


def summarize_gene_all_to_all(gene_payload):
    """Per-gene worker (module-level so ProcessPoolExecutor can pickle it)."""
    gene_id = gene_payload["geneID"]
    group1_bias = gene_payload["group1_bias"]
    group2_bias = gene_payload["group2_bias"]
    group1 = gene_payload["group1"]
    group2 = gene_payload["group2"]
    rows = gene_payload["rows"]

    best_all = {event_type: None for event_type in EVENT_TYPES}
    best_between = {event_type: None for event_type in EVENT_TYPES}

    for tx_a, tx_b in combinations(rows, 2):
        event_flags = pairwise_as_detection(tx_a, tx_b)
        # Pair inherits the less-annotated member's tier
        pair_tier = max(tx_a["anno_tier"], tx_b["anno_tier"])
        is_between_group_pair = (
            (tx_a["sexClass_combined"] == group1_bias and tx_b["sexClass_combined"] == group2_bias)
            or (tx_a["sexClass_combined"] == group2_bias and tx_b["sexClass_combined"] == group1_bias)
        )
        for event_type, event_present in event_flags.items():
            if not event_present:
                continue
            if best_all[event_type] is None or pair_tier < best_all[event_type]:
                best_all[event_type] = pair_tier
            if is_between_group_pair:
                if best_between[event_type] is None or pair_tier < best_between[event_type]:
                    best_between[event_type] = pair_tier

    result_row = {"geneID": gene_id}
    for event_type in EVENT_TYPES:
        output_name = EVENT_OUTPUT[event_type]
        result_row[output_name] = TIER_OUTPUT.get(best_all[event_type], "none")
        result_row[f"{output_name}_{group1}v{group2}"] = TIER_OUTPUT.get(best_between[event_type], "none")
    return result_row


def main():
    parser = argparse.ArgumentParser(description="All-to-all AS detection collapsed to gene-level AS columns.")
    parser.add_argument("-i", "--input", dest="ujc_file", required=True)
    parser.add_argument("-e", "--erp-input", dest="erp_file", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-g1", "--group1", required=True)
    parser.add_argument("-g2", "--group2", required=True)
    parser.add_argument("--group-class-col", default="sexClass")
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()

    cpus = max(1, args.cpus)
    group1_bias = f"{args.group1}_bias"
    group2_bias = f"{args.group2}_bias"

    ujc_df = pd.read_csv(args.ujc_file, low_memory=False)
    ujc_df = ujc_df.rename(columns={args.group_class_col: "sexClass"})

    # Merge ERP_plus groupClass as a fallback bias label per UJC.
    erp_df = pd.read_csv(args.erp_file, low_memory=False)
    erp_group_class = (
        erp_df.groupby(["geneID", "ERP_plus"], sort=False)["sexClass"]
        .agg(lambda s: group1_bias if (s == group1_bias).any()
             else group2_bias if (s == group2_bias).any()
             else next((str(v) for v in s if pd.notna(v) and str(v) != ""), ""))
        .reset_index(name="sexClass_ERPp")
    )
    ujc_df = pd.merge(ujc_df, erp_group_class, on=["geneID", "ERP_plus"], how="left")

    ujc_class = ujc_df["sexClass"].fillna("")
    erp_class = ujc_df["sexClass_ERPp"].fillna("")
    ujc_df["sexClass_combined"] = ujc_class.where(
        ujc_class.isin([group1_bias, group2_bias]),
        erp_class.where(erp_class.isin([group1_bias, group2_bias]), ujc_class),
    )

    # Eligible = non-ISM (has a tier) and multi-exon-region
    ujc_df["anno_tier"] = ujc_df["annotation_status"].map(ANNOTATION_TIER_MAP)
    ujc_df["num_exon_regions"] = ujc_df["ERP"].map(lambda e: len(erp_to_set(str(e))))
    eligible_mask = ujc_df["anno_tier"].notna() & (ujc_df["num_exon_regions"] >= 2)

    needed = ["geneID", "jxnHash", "ERP", "IR_ER", "dataOnlyER_ID", "anno_tier", "sexClass_combined"]
    eligible_df = ujc_df.loc[eligible_mask, needed].copy()
    eligible_df["IR_ER"] = eligible_df["IR_ER"].fillna("")
    eligible_df["dataOnlyER_ID"] = eligible_df["dataOnlyER_ID"].fillna("")
    eligible_df["anno_tier"] = eligible_df["anno_tier"].astype(int)

    payloads = []
    for gene_id, gene_rows in eligible_df.groupby("geneID", sort=False):
        payloads.append({
            "geneID": gene_id,
            "rows": gene_rows[needed[1:]].to_dict(orient="records"),
            "group1": args.group1,
            "group2": args.group2,
            "group1_bias": group1_bias,
            "group2_bias": group2_bias,
        })

    print(f"Genes with eligible UJCs: {len(payloads):,}", flush=True)
    if cpus <= 1:
        as_results = [summarize_gene_all_to_all(p) for p in payloads]
    else:
        with ProcessPoolExecutor(max_workers=cpus) as executor:
            as_results = list(executor.map(summarize_gene_all_to_all, payloads))

    # Genes with no eligible pair get 'none' for every event column.
    none_row = {}
    for event_type in EVENT_TYPES:
        output_name = EVENT_OUTPUT[event_type]
        none_row[output_name] = "none"
        none_row[f"{output_name}_{args.group1}v{args.group2}"] = "none"
    seen = {row["geneID"] for row in as_results}
    for gene_id in ujc_df["geneID"].drop_duplicates():
        if gene_id not in seen:
            as_results.append({"geneID": gene_id, **none_row})

    as_df = pd.DataFrame(as_results)
    as_df.to_csv(args.output, index=False)

    print(f"Species: {args.species}", flush=True)
    print(f"Genes in output: {as_df.shape[0]:,}", flush=True)
    print(f"Output saved to: {args.output}", flush=True)


if __name__ == "__main__":
    main()
