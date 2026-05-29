#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 16:48:56 2025

@author: nkeil
"""
import pandas as pd
import os
import matplotlib.pyplot as plt

# Load the datasets 
ind = '/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts'

cd4_df = pd.read_csv(os.path.join(ind,"flag_analyze_frag_DE_CD4.csv"), sep=",")  # CD4 dataset
cd8_df = pd.read_csv(os.path.join(ind,"flag_analyze_frag_DE_CD8.csv"), sep=",")  # CD8 dataset

# Merge datasets on 'featureID' to compare 'flag_analyze_frag_DE' values
merged_df = cd4_df[['featureID', 'flag_analyze_frag_DE']].merge(
    cd8_df[['featureID', 'flag_analyze_frag_DE']],
    on='featureID',
    suffixes=('_CD4', '_CD8')
)

# Compute flag_differentially_detected (1 if different, 0 if same)
merged_df['flag_differentially_detected'] = (
    merged_df['flag_analyze_frag_DE_CD4'] != merged_df['flag_analyze_frag_DE_CD8']
).astype(int)


# Extract gene ID (everything before the first colon in featureID)
merged_df['geneID'] = merged_df['featureID'].apply(lambda x: x.split(":")[0])

#Load t1d and immune gene lists
lst_dir = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists"
t1d_genes = set(pd.read_csv(os.path.join(lst_dir,"T1D_candidate_genes_robertson_2021_ensembl.txt"))["ensembl_geneID"])
immune_genes = set(pd.read_csv(os.path.join(lst_dir,"immunobase_hg19_ensembl_entrez.csv"))["Ensembl_ID"])

### **Function to count analyzable fragments per gene category**
def count_analyzable(df, gene_list=None):
    if gene_list is not None:
        df = df[df['geneID'].isin(gene_list)]  # Filter by gene list

    # Count fragments that are analyzable in different conditions
    cd4_analyzable = df[df['flag_analyze_frag_DE_CD4'] == 1].shape[0]
    cd8_analyzable = df[df['flag_analyze_frag_DE_CD8'] == 1].shape[0]
    both_analyzable = df[(df['flag_analyze_frag_DE_CD4'] == 1) & (df['flag_analyze_frag_DE_CD8'] == 1)].shape[0]
    cd4_only_analyzable = df[(df['flag_analyze_frag_DE_CD4'] == 1) & (df['flag_analyze_frag_DE_CD8'] == 0)].shape[0]
    cd8_only_analyzable = df[(df['flag_analyze_frag_DE_CD4'] == 0) & (df['flag_analyze_frag_DE_CD8'] == 1)].shape[0]

    return [cd4_analyzable, cd8_analyzable, both_analyzable, cd4_only_analyzable, cd8_only_analyzable]

# Compute analyzable fragment counts for all genes, immune genes, and T1D genes
counts_all = count_analyzable(merged_df)
counts_immune = count_analyzable(merged_df, immune_genes)
counts_t1d = count_analyzable(merged_df, t1d_genes)

# Create a summary table
analyzable_summary_table = pd.DataFrame({
    "Category": ["CD4 analyzable", "CD8 analyzable", "Analyzable in both", "CD4 only analyzable", "CD8 only analyzable"],
    "All Genes": counts_all,
    "Immune Genes": counts_immune,
    "T1D Genes": counts_t1d
})

analyzable_summary_table.to_csv(os.path.join(ind,'analyzable_frag_counts.csv'), index = False)


# Identify genes with at least one differentially detected fragment
diff_genes = set(merged_df[merged_df['flag_differentially_detected'] == 1]['geneID'])

# Count differentially detected genes for each category
num_diff_all = len(diff_genes)
num_diff_t1d = len(diff_genes.intersection(t1d_genes))
num_diff_immune = len(diff_genes.intersection(immune_genes))

#!!! Total genes in each category
# total_genes = len(set(merged_df['geneID']))
# total_t1d = len(t1d_genes)
# total_immune = len(immune_genes)

total_genes = 6994
total_t1d = 61
total_immune = 815

# Compute proportions
prop_diff_all = num_diff_all / total_genes if total_genes > 0 else 0
prop_diff_t1d = num_diff_t1d / total_t1d if total_t1d > 0 else 0
prop_diff_immune = num_diff_immune / total_immune if total_immune > 0 else 0

# Create the gene-level summary table
gene_summary_table = pd.DataFrame({
    "category": ["All Genes", "Immune Genes", "T1D Genes"],
    "total_genes": [total_genes, total_immune, total_t1d],
    "num_genes_w_diff_det_frags": [num_diff_all, num_diff_immune, num_diff_t1d],
    "prop_genes_w_diff_det_frags": [prop_diff_all, prop_diff_immune, prop_diff_t1d]
})

gene_summary_table.to_csv(os.path.join(ind,'count_genes_w_diff_det_frags.csv'), index = False)



## Function to generate and save pie charts
#def create_pie_chart(detected, total, title, filename):
#     plt.figure(figsize=(6, 6))
#     labels = ['Differentially Detected', 'Not Differentially Detected']
#     sizes = [detected, total - detected]
#     plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
#     plt.title(title)
#     plt.savefig(filename)
#     plt.show()

# # Generate pie charts
# create_pie_chart(num_diff_all, len(set(merged_df['geneID'])), "All Genes Differentially Detected", "all_genes.pdf")
# create_pie_chart(num_diff_immune, total_immune, "Immune Genes Differentially Detected", "immune_genes.pdf")
# create_pie_chart(num_diff_t1d, total_t1d, "T1D Genes Differentially Detected", "t1d_genes.pdf")

import seaborn as sns

# Match category order used in your DE-genes plot so colors line up 1:1
cat_order = ["All Genes", "Immune Genes", "T1D Genes"]

# Build the same Set2 palette and bind it to categories
palette = sns.color_palette("Set2", n_colors=len(cat_order))
cat2color = dict(zip(cat_order, palette))

# Prepare fragment proportions from your existing gene_summary_table
frag_df = gene_summary_table[["category", "prop_genes_w_diff_det_frags"]].copy()
frag_df.columns = ["Category", "Proportion"]
frag_df["Category"] = pd.Categorical(frag_df["Category"], categories=cat_order, ordered=True)
frag_df = frag_df.sort_values("Category")

plt.figure(figsize=(6, 4))
ax = sns.barplot(
    data=frag_df,
    x="Category",
    y="Proportion",
    palette=[cat2color[c] for c in frag_df["Category"]]
)

ax.set_ylim(0, 1)
ax.set_ylabel("Proportion of Genes with\nDifferentially Detected Fragments")
ax.set_xlabel("Gene Category")
plt.tight_layout()

# Save without any bar labels
plt.savefig(os.path.join(ind, "prop_genes_diff_detected_barplot.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(ind, "prop_genes_diff_detected_barplot.png"), dpi=300, bbox_inches="tight")
plt.close()

