# -*- coding: utf-8 -*-


import os
import pandas as pd

# Set paths
ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"

# 1) Fragment flags file (the one with CD4/CD8 fragment info)
frag_flag_file = os.path.join(ind, "quantify_t1d_pacbio_transcripts/flag_analyzable_CD4_CD8_combined.csv")  

# 2) Intron annotation file (with ef_ir_flag)
intron_annot_file = os.path.join(ind,"allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv" )

# 3) Gene lists (same as in your existing script)
lst_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists"
t1d_list_file = os.path.join(lst_dir, "T1D_candidate_genes_robertson_2021_ensembl.txt")
immune_list_file = os.path.join(lst_dir, "immunobase_hg19_ensembl_entrez.csv")

# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------

# Fragment-level table with CD4/CD8 flags
frag_df = pd.read_csv(frag_flag_file, sep=",")

# Expect at least these columns:
#   featureID, geneID,
#   flag_frag_num_on_cd4, flag_frag_num_on_cd8
required_cols = [
    "featureID",
    "geneID",
    "flag_analyze_frag_DE_cd4",
    "flag_analyze_frag_DE_cd8",
]
missing = [c for c in required_cols if c not in frag_df.columns]
if missing:
    raise ValueError(f"Missing expected columns in fragment file: {missing}")

# Intron annotation table (ef_ir_flag)
intron_df = pd.read_csv(intron_annot_file, sep=",")

# Expect columns: ef_id, ef_ir_flag (plus others)
if "ef_id" not in intron_df.columns or "ef_ir_flag" not in intron_df.columns:
    raise ValueError("Intron file must contain 'ef_id' and 'ef_ir_flag' columns.")

# Keep only EFs that are intron-retention (ef_ir_flag == 1)
ir_introns = intron_df[intron_df["ef_ir_flag"] == 1].copy()
ir_ef_ids = set(ir_introns["ef_id"])

# ---------------------------------------------------------------------
# Restrict fragment table to IR fragments only
# ---------------------------------------------------------------------
frag_ir_df = frag_df[frag_df["featureID"].isin(ir_ef_ids)].copy()

# Define whether each IR EF is "on" in CD4/CD8
frag_ir_df["cd4_analyzable"] = frag_ir_df["flag_analyze_frag_DE_cd4"] == 1
frag_ir_df["cd8_analyzable"] = frag_ir_df["flag_analyze_frag_DE_cd8"] == 1
frag_ir_df["analyzable_either"] = frag_ir_df["cd4_analyzable"] | frag_ir_df["cd8_analyzable"]

# ---------------------------------------------------------------------
# Identify IR genes (at least one IR EF is "on" in CD4 or CD8)
# ---------------------------------------------------------------------
ir_gene_ids = set(frag_ir_df.loc[frag_ir_df["analyzable_either"], "geneID"])

# ---------------------------------------------------------------------
# Candidate gene lists
# ---------------------------------------------------------------------
t1d_genes = set(pd.read_csv(t1d_list_file)["ensembl_geneID"])
immune_genes = set(pd.read_csv(immune_list_file)["Ensembl_ID"])

# ---------------------------------------------------------------------
# Counts of IR genes in each category
# ---------------------------------------------------------------------
num_ir_all = len(ir_gene_ids)
num_ir_immune = len(ir_gene_ids.intersection(immune_genes))
num_ir_t1d = len(ir_gene_ids.intersection(t1d_genes))

# ---------------------------------------------------------------------
# Total genes in each category (use same denominators as your DE script)
# ---------------------------------------------------------------------
# If you want to recompute these automatically instead, comment these
# lines out and uncomment the alternative block below.
total_genes = 6994
total_t1d = 61
total_immune = 815

# Alternative (automatic) totals if desired:
# total_genes = len(set(frag_df["geneID"]))
# total_t1d = len(t1d_genes)
# total_immune = len(immune_genes)

# Proportions
prop_ir_all = num_ir_all / total_genes if total_genes > 0 else 0
prop_ir_immune = num_ir_immune / total_immune if total_immune > 0 else 0
prop_ir_t1d = num_ir_t1d / total_t1d if total_t1d > 0 else 0

# ---------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------
ir_gene_summary_table = pd.DataFrame({
    "category": ["All Genes", "Immune Genes", "T1D Genes"],
    "total_genes": [total_genes, total_immune, total_t1d],
    "num_genes_w_IR_frags_on": [num_ir_all, num_ir_immune, num_ir_t1d],
    "prop_genes_w_IR_frags_on": [prop_ir_all, prop_ir_immune, prop_ir_t1d]
})

# Output
out_file = os.path.join(ind, "count_genes_w_IR_frags_on.csv")
ir_gene_summary_table.to_csv(out_file, index=False)
