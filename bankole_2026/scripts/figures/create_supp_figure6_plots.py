#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 16:57:09 2026

@author: mgaran
"""

import pandas as pd
import matplotlib.pyplot as plt

# Project paths.
PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
DUP = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/duplication_vs_splicing/analysis_adp"
TABLES = f"{PROJ}/Tables"
FIGURES = f"{PROJ}/Figures"


# =============================================================================
# D. melanogaster
# =============================================================================

dup_dmel = pd.read_csv(f"{DUP}/dmel_singleton_vs_dup.csv")
data_dmel = pd.read_csv(
    f"{TABLES}/gene_level_summary_from_dataFile_jxnHash_w_flags_dmel.csv"
)

# Number of non-ISM transcripts per gene.
data_dmel["num_nonISM_UJC"] = (
    data_dmel["num_UJCanno_UJC"]
    + data_dmel["num_ERPanno_UJC"]
    + data_dmel["num_ERPnovel_UJC"]
)

# Match the copy-number gene ID column to the gene-summary gene ID column.
dup_dmel = dup_dmel.rename(columns={"dmel": "geneID"})

merged_dmel = pd.merge(
    dup_dmel,
    data_dmel,
    on="geneID",
    how="outer",
    indicator="merge_check",
)

print("dmel merge_check:")
print(merged_dmel["merge_check"].value_counts())

# Keep genes present in both files.
merged_dmel = merged_dmel[merged_dmel["merge_check"] == "both"]

x_dmel = merged_dmel["num_FBGNs"].to_numpy(dtype=float)
y_dmel = merged_dmel["num_nonISM_UJC"].to_numpy(dtype=float)
y2_dmel = merged_dmel["num_UJCanno_UJC"].to_numpy(dtype=float)

# Plot all non-ISM transcripts.
plt.figure()
plt.scatter(x_dmel, y_dmel, s=30, alpha=0.6, color="#D4B08C")
plt.title("dmel copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/copy_number_vs_num_nonISM_UJC_scatter_dmel.png",
    bbox_inches="tight",
)
plt.close()

# Plot annotated transcripts.
plt.figure()
plt.scatter(x_dmel, y2_dmel, s=30, alpha=0.6, color="#966729")
plt.title("dmel copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of annotated transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/supp_figure6dmel.png",
    bbox_inches="tight",
)
plt.close()


# =============================================================================
# D. simulans
# =============================================================================

dup_dsim = pd.read_csv(f"{DUP}/dsim_5species_Fbgn_2_Hahn_GLEANR.csv")
data_dsim = pd.read_csv(
    f"{TABLES}/gene_level_summary_from_dataFile_jxnHash_w_flags_dsim.csv"
)

# Number of non-ISM transcripts per gene.
data_dsim["num_nonISM_UJC"] = (
    data_dsim["num_UJCanno_UJC"]
    + data_dsim["num_ERPanno_UJC"]
    + data_dsim["num_ERPnovel_UJC"]
)

# Match the gene-summary gene ID column to the copy-number file.
data_dsim = data_dsim.rename(columns={"geneID": "dsim_FBgn"})

merged_dsim = pd.merge(
    dup_dsim,
    data_dsim,
    on="dsim_FBgn",
    how="outer",
    indicator="merge_check",
)

print("dsim merge_check:")
print(merged_dsim["merge_check"].value_counts())

# Keep genes present in both files.
merged_dsim = merged_dsim[merged_dsim["merge_check"] == "both"]

x_dsim = merged_dsim["num_GLEANRIDs_sim"].to_numpy(dtype=float)
y_dsim = merged_dsim["num_nonISM_UJC"].to_numpy(dtype=float)
y2_dsim = merged_dsim["num_UJCanno_UJC"].to_numpy(dtype=float)

# Plot all non-ISM transcripts.
plt.figure()
plt.scatter(x_dsim, y_dsim, s=30, alpha=0.6, color="#9FC1E8")
plt.title("dsim copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/copy_number_vs_num_nonISM_UJC_scatter_dsim.png",
    bbox_inches="tight",
)
plt.close()

# Plot annotated transcripts.
plt.figure()
plt.scatter(x_dsim, y2_dsim, s=30, alpha=0.6, color="#3F78C1")
plt.title("dsim copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of annotated transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/supp_figure6dsim.png",
    bbox_inches="tight",
)
plt.close()


# =============================================================================
# D. yakuba
# =============================================================================

dup_dyak = pd.read_csv(f"{DUP}/dyak_Hahn_NCBI_GLEANR_link.csv")
data_dyak = pd.read_csv(
    f"{TABLES}/gene_level_summary_from_dataFile_jxnHash_w_flags_dyak.csv"
)

# Number of non-ISM transcripts per gene.
data_dyak["num_nonISM_UJC"] = (
    data_dyak["num_UJCanno_UJC"]
    + data_dyak["num_ERPanno_UJC"]
    + data_dyak["num_ERPnovel_UJC"]
)

# Convert the NCBI gene IDs to the LOC IDs used in the gene-summary file.
dup_dyak = dup_dyak.rename(columns={"NCBI_geneid": "geneID"})
dup_dyak["geneID"] = dup_dyak["geneID"].apply(lambda x: str(int(x)))
dup_dyak["geneID"] = "LOC" + dup_dyak["geneID"].astype(str)

merged_dyak = pd.merge(
    dup_dyak,
    data_dyak,
    on="geneID",
    how="outer",
    indicator="merge_check",
)

print("dyak merge_check:")
print(merged_dyak["merge_check"].value_counts())

# Keep genes present in both files.
merged_dyak = merged_dyak[merged_dyak["merge_check"] == "both"]

x_dyak = merged_dyak["num_GLEANRIDs_yak"].to_numpy(dtype=float)
y_dyak = merged_dyak["num_nonISM_UJC"].to_numpy(dtype=float)
y2_dyak = merged_dyak["num_UJCanno_UJC"].to_numpy(dtype=float)

# Plot all non-ISM transcripts.
plt.figure()
plt.scatter(x_dyak, y_dyak, s=30, alpha=0.6, color="#B8B9BA")
plt.title("dyak copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/copy_number_vs_num_nonISM_UJC_scatter_dyak.png",
    bbox_inches="tight",
)
plt.close()

# Plot annotated transcripts.
plt.figure()
plt.scatter(x_dyak, y2_dyak, s=30, alpha=0.6, color="#717273")
plt.title("dyak copy number vs number transcripts")
plt.xlabel("copy number")
plt.ylabel("number of annotated transcripts")
plt.xlim(0, 50)
plt.savefig(
    f"{FIGURES}/supp_figure6dyak.png",
    bbox_inches="tight",
)
plt.close()

print("Supplementary Figure 6 plots complete.")