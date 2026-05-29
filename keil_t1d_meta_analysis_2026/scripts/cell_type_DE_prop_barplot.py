#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 21:21:50 2025

@author: nkeil
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Set output directory for pdf figures
outd = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/manuscript/pdf_figs"

#Plot proportion of genes with sig cell type difference

#20% of frags in gene significant with p value of 0.001

ind = "/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_controls_CD4_CD8/DE_gene_counts"

infile = "gene_DE_cnts_cellType_effect_sig_0.001_cut_0.2.csv"

inpath = os.path.join(ind,infile)

df = pd.read_csv(inpath)

# Plot with different colors for each bar
palette = sns.color_palette("Set2", n_colors=3)

plt.figure(figsize=(6, 4))
sns.barplot(data=df, x='Category', y='Prop_genes_DE', palette=palette)

plt.ylabel("Proportion of Genes Differentially Expressed")
plt.xlabel("Gene Category")
#plt.title("Proportion of DE Genes by Category")
plt.ylim(0, 1)
plt.tight_layout()
#plt.show()

# Save as PDF and PNG (before plt.show or plt.close!)
plt.savefig(f"{outd}/prop_genes_de_plot.pdf", format="pdf", bbox_inches="tight")
plt.savefig(f"{outd}/prop_genes_de_plot.png", format="png", bbox_inches="tight")

plt.close()  