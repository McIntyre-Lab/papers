#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Counts features significant for flag_case OR (femal*case interaction) and reports
counts/proportions for:
  - all genes
  - immune genes
  - T1D candidate genes
for each cell in ["CD4","CD8"].

Also reports gene-level tested and significant counts/proportions:
  - num_genes_tested
  - num_genes_with_sig_feature
  - prop_genes_with_sig_feature

- Performs counts separately for EXONS and INTRONS using an annotation file
  containing ef_id and ef_ir_flag (0 exon, 1 intron).
- Writes two separate output files.
"""

import os
import pandas as pd


# -------------------- USER SETTINGS --------------------
ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_age_sex"
sig = 0.05

list_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists"

immune_genes_file = os.path.join(list_dir, "immunobase_hg19_ensembl_entrez.csv")
immune_genes_col  = "Ensembl_ID"
immune_sep = ","  # change to "\t" if needed

t1d_genes_file = os.path.join(list_dir, "T1D_candidate_genes_robertson_2021_ensembl.txt")
t1d_genes_col  = "ensembl_geneID"   # set to the real column name in that file
t1d_sep = "\t"

# Annotation file mapping each fragment (ef_id) to ef_ir_flag (0 exon / 1 intron)
# MUST contain columns: ef_id, ef_ir_flag
frag_annot_file = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv"  
frag_annot_sep = ","  # change to "\t" if needed

out_csv_exons  = os.path.join(ind, "frags_case_effect_counts_EXONS.csv")
out_csv_introns = os.path.join(ind, "frags_case_effect_counts_INTRONS.csv")
# -------------------------------------------------------


def read_ensembl_list_table(path, col, sep):
    """Read Ensembl IDs from a delimited table (CSV/TSV) column."""
    df = pd.read_csv(path, sep=sep)
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in {path}. Columns: {list(df.columns)}")
    ids = (
        df[col]
        .astype(str)
        .str.strip()
        .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
        .dropna()
        .tolist()
    )
    return set(ids)


def clean_probf(series):
    """
    Convert ProbF to numeric, handling values like "<.0000..." by turning them into
    "0.0000..." (replace leading '<' with '0').
    Unparseable values become NaN.
    """
    s = series.astype(str).str.strip()
    s = s.str.replace(r"^\<\s*", "0", regex=True)  # "<.0001" -> "0.0001"
    return pd.to_numeric(s, errors="coerce")


def load_fragment_annotation(path, sep):
    """
    Load fragment annotation mapping ef_id -> ef_ir_flag (0 exon, 1 intron).
    Returns a df with columns: featureID, ef_ir_flag.
    """
    ann = pd.read_csv(path, sep=sep)

    needed = {"ef_id", "ef_ir_flag"}
    missing = needed - set(ann.columns)
    if missing:
        raise ValueError(f"Annotation file missing columns {missing}. Found columns: {list(ann.columns)}")

    ann = ann[["ef_id", "ef_ir_flag"]].copy()
    ann["ef_id"] = ann["ef_id"].astype(str).str.strip()
    ann["ef_ir_flag"] = pd.to_numeric(ann["ef_ir_flag"], errors="coerce")

    # ef_id should correspond to featureID in DE table
    ann = ann.rename(columns={"ef_id": "featureID"})

    # If duplicates exist, keep first (they should be identical; warn if not)
    if ann["featureID"].duplicated().any():
        # optional: check if conflicting flags exist
        conflict = (
            ann.groupby("featureID")["ef_ir_flag"].nunique(dropna=False)
              .reset_index(name="n_unique_flags")
              .query("n_unique_flags > 1")
        )
        if len(conflict) > 0:
            print(f"WARNING: {len(conflict)} featureIDs have conflicting ef_ir_flag values in annotation file.")
        ann = ann.drop_duplicates("featureID", keep="first")

    return ann


def compute_counts(df_sub, cell, gene_set_name, sig, subset_label):
    """
    df_sub must already be filtered to:
      - effects of interest
      - a particular gene set (all/immune/t1d)
      - and exon/intron subset

    Returns one row dict with both feature- and gene-level metrics.
    """
    # Significant rows (ignore NaN ProbF)
    df_sub["is_sig"] = df_sub["ProbF"].notna() & (df_sub["ProbF"] < sig)

    # -------- Feature-level metrics --------
    feature_sig = df_sub.groupby("featureID")["is_sig"].any()
    n_features_sig = int(feature_sig.sum())
    n_features_tested = int(feature_sig.shape[0])
    prop_features_sig = n_features_sig / n_features_tested if n_features_tested else float("nan")

    # -------- Gene-level metrics --------
    num_genes_tested = int(df_sub["geneID"].nunique()) if len(df_sub) else 0

    if n_features_tested:
        feat_to_gene = df_sub.drop_duplicates("featureID").set_index("featureID")["geneID"]
        sig_features = feature_sig[feature_sig].index
        num_genes_with_sig_feature = int(feat_to_gene.loc[sig_features].nunique()) if len(sig_features) else 0
    else:
        num_genes_with_sig_feature = 0

    prop_genes_with_sig_feature = (
        num_genes_with_sig_feature / num_genes_tested if num_genes_tested else float("nan")
    )

    return {
        "cell": cell,
        "subset": subset_label,  # EXONS / INTRONS
        "gene_set": gene_set_name,
        "alpha": sig,

        "num_features_tested": n_features_tested,
        "num_features_sig_case_or_interaction": n_features_sig,
        "prop_features_sig_case_or_interaction": prop_features_sig,

        "num_genes_tested": num_genes_tested,
        "num_genes_with_sig_feature": num_genes_with_sig_feature,
        "prop_genes_with_sig_feature": prop_genes_with_sig_feature,
    }


def main():
    immune_genes = read_ensembl_list_table(immune_genes_file, col=immune_genes_col, sep=immune_sep)
    t1d_genes = read_ensembl_list_table(t1d_genes_file, col=t1d_genes_col, sep=t1d_sep)

    ann = load_fragment_annotation(frag_annot_file, sep=frag_annot_sep)

    rows_exons = []
    rows_introns = []

    for cell in ["CD4", "CD8"]:
        infile = os.path.join(ind, f"de_t1_frag_{cell}_FbyCase.csv")
        df = pd.read_csv(infile)

        # Clean ProbF FIRST so subsets inherit numeric dtype
        df["ProbF"] = clean_probf(df["ProbF"])

        # Clean Effect for robust matching
        df["Effect_clean"] = df["Effect"].astype(str).str.strip().str.lower()

        # Extract geneID from featureID like "ENSG...:ER1:EF1"
        df["geneID"] = df["featureID"].astype(str).str.split(":", n=1).str[0]

        # Merge in ef_ir_flag (0 exon, 1 intron)
        df = df.merge(ann, on="featureID", how="left")

        n_missing_flag = int(df["ef_ir_flag"].isna().sum())
        if n_missing_flag > 0:
            print(f"WARNING ({cell}): {n_missing_flag} rows missing ef_ir_flag after merge. They will be dropped for exon/intron summaries.")

        # Keep only rows that have a valid exon/intron flag
        df = df[df["ef_ir_flag"].isin([0, 1])].copy()

        # Effects of interest (include both spellings if they exist in your files)
        effects_of_interest = {"flag_case", "flag_femal*flag_case", "flag_female*flag_case"}

        # Keep only these effects
        sub0 = df[df["Effect_clean"].isin(effects_of_interest)].copy()

        gene_sets = [
            ("all_genes", None),
            ("immune_genes", immune_genes),
            ("t1d_genes", t1d_genes),
        ]

        # Split exon/intron
        sub_exons = sub0[sub0["ef_ir_flag"] == 0].copy()
        sub_introns = sub0[sub0["ef_ir_flag"] == 1].copy()

        for gene_set_name, gene_set in gene_sets:
            ex = sub_exons if gene_set is None else sub_exons[sub_exons["geneID"].isin(gene_set)].copy()
            intr = sub_introns if gene_set is None else sub_introns[sub_introns["geneID"].isin(gene_set)].copy()

            rows_exons.append(compute_counts(ex, cell, gene_set_name, sig, subset_label="EXONS"))
            rows_introns.append(compute_counts(intr, cell, gene_set_name, sig, subset_label="INTRONS"))

    df_exons = pd.DataFrame(rows_exons)
    df_introns = pd.DataFrame(rows_introns)

    print("\n=== EXONS ===")
    print(df_exons.to_string(index=False))
    df_exons.to_csv(out_csv_exons, index=False)
    print(f"\nWrote: {out_csv_exons}")

    print("\n=== INTRONS ===")
    print(df_introns.to_string(index=False))
    df_introns.to_csv(out_csv_introns, index=False)
    print(f"\nWrote: {out_csv_introns}")


if __name__ == "__main__":
    main()
