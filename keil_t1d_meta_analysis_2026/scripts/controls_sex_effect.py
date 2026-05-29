#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 13:08:27 2023

@author: nkeil
"""

import pandas as pd
import numpy as np
import argparse

inDir="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/quantify_t1d_pacbio_transcripts/model_DE_frag_controls_CD4_CD8/"

t1_path=inDir+"de_t1_frag_control.csv"

t1DF=pd.read_csv(t1_path, sep=',', header=0, low_memory=False)
t1DF['geneID'] = t1DF['featureID'].str.split(':').str[0]

#Dealing with <0.0000000001

t1DF['ProbF'] = pd.to_numeric(t1DF['ProbF'].str.replace('<', ''), errors='coerce')

# Creating 'flag_sig_0.001' column: 1 if 'ProbF' < 0.001, else 0
t1DF['flag_sig_001'] = (t1DF['ProbF'] < 0.001).astype(int)



# Creating 'flag_sig_0.00001' column: 1 if 'ProbF' < 0.00001, else 0
t1DF['flag_sig_00001'] = (t1DF['ProbF'] < 0.00001).astype(int)

##Sex effect



# Step 1: Group by both 'geneID' and 'Effect' and aggregate the data
gene_effect_DF = t1DF.groupby(['geneID', 'Effect']).agg(
    cnt_frags_sig_001=('flag_sig_001', 'sum'),  # Sum of frags significant at 0.001
    cnt_frags_sig_00001=('flag_sig_00001', 'sum'),  # Sum of frags significant at 0.00001
    total_frags=('geneID', 'size')  # Count of total frags
)

# Step 2: Calculate the proportions
gene_effect_DF['prop_frags_sig_001'] = gene_effect_DF['cnt_frags_sig_001'] / gene_effect_DF['total_frags']
gene_effect_DF['prop_frags_sig_00001'] = gene_effect_DF['cnt_frags_sig_00001'] / gene_effect_DF['total_frags']

# Reset index to make 'geneID' and 'Effect' columns again
gene_effect_DF.reset_index(inplace=True)

gene_sex_effect_DF=gene_effect_DF[gene_effect_DF['Effect'] == 'flag_female']

import seaborn as sns
import matplotlib.pyplot as plt

filter_prop_001_DF=gene_sex_effect_DF[gene_sex_effect_DF['prop_frags_sig_001']> 0]

filter_prop_00001_DF=gene_sex_effect_DF[gene_sex_effect_DF['prop_frags_sig_00001']> 0]

#Hist - no.of genes with rpoportion of genes significant

# Setting up the plot dimensions and style
plt.figure(figsize=(12, 6))
plt.style.use('ggplot')

# Histogram for proportion_of_frags_sig_0_001
plt.hist(filter_prop_001_DF['prop_frags_sig_001'], bins=20, alpha=0.5, label='Proportion SIG 0.001')

# Histogram for proportion_of_frags_sig_0_00001
plt.hist(filter_prop_00001_DF['prop_frags_sig_00001'], bins=20, alpha=0.5, label='Proportion SIG 0.00001')

# Adding labels and title
plt.title('Distribution of Proportions - Sex Effect')
plt.xlabel('Proportion Value')
plt.ylabel('Number of Genes')

# Adding a legend
plt.legend()

#Save the figure
plt.savefig(inDir+"hist_prop_sex_effect_frags_DE.png")
plt.show()
plt.close()


#From plots we use cut off of 0.2 as marker for DE

# Flag for proportion_of_frags_sig_0_001 greater than 0.2
gene_sex_effect_DF['flag_gene_DE_001'] = (gene_sex_effect_DF['prop_frags_sig_001'] > 0.2).astype(int)

# Flag for proportion_of_frags_sig_0_00001 greater than 0.2
gene_sex_effect_DF['flag_gene_DE_00001'] = (gene_sex_effect_DF['prop_frags_sig_00001'] > 0.2).astype(int)

#Import immune genes and T1D gene lists

t1d_gene="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"
immune_gene="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists/immunobase_hg19_ensembl_entrez.csv"

t1d_geneDF=pd.read_csv(t1d_gene, sep=',', header=0, low_memory=False)
immune_geneDF=pd.read_csv(immune_gene, sep=',', header=0, low_memory=False)

#Flag immune genes and t1d genes

gene_sex_effect_DF['flag_T1D_gene'] = gene_sex_effect_DF['geneID'].isin(t1d_geneDF['ensembl_geneID']).astype(int)
gene_sex_effect_DF['flag_immune_gene'] = gene_sex_effect_DF['geneID'].isin(immune_geneDF['Ensembl_ID']).astype(int)

#Cnt DE genes
cnt_DE_all_001=gene_sex_effect_DF['flag_gene_DE_001'].sum()
cnt_DE_all_00001=gene_sex_effect_DF['flag_gene_DE_00001'].sum()
cnt_DE_immune_001=(gene_sex_effect_DF['flag_gene_DE_001'] * gene_sex_effect_DF['flag_immune_gene']).sum()
cnt_DE_immune_0001=(gene_sex_effect_DF['flag_gene_DE_00001'] * gene_sex_effect_DF['flag_immune_gene']).sum()
cnt_DE_t1d_001=(gene_sex_effect_DF['flag_gene_DE_001'] * gene_sex_effect_DF['flag_T1D_gene']).sum()
cnt_DE_t1d_00001=(gene_sex_effect_DF['flag_gene_DE_00001'] * gene_sex_effect_DF['flag_T1D_gene']).sum()

#Get counts of analyzable genes
cnt_all=len(gene_sex_effect_DF)
cnt_immune=gene_sex_effect_DF['flag_immune_gene'].sum()
cnt_T1D=gene_sex_effect_DF['flag_T1D_gene'].sum()

# Calculate proportions for 0.001
prop_DE_all_001 = cnt_DE_all_001 / cnt_all
prop_DE_immune_001 = cnt_DE_immune_001 / cnt_immune if cnt_immune > 0 else 0
prop_DE_t1d_001 = cnt_DE_t1d_001 / cnt_T1D if cnt_T1D > 0 else 0

# Calculate proportions for 0.00001
prop_DE_all_00001 = cnt_DE_all_00001 / cnt_all
prop_DE_immune_00001 = cnt_DE_immune_0001 / cnt_immune if cnt_immune > 0 else 0
prop_DE_t1d_00001 = cnt_DE_t1d_00001 / cnt_T1D if cnt_T1D > 0 else 0

# Proportions data for plotting with the order of Immune and T1D genes switched
proportions_001 = [prop_DE_all_001, prop_DE_immune_001, prop_DE_t1d_001]
proportions_00001 = [prop_DE_all_00001, prop_DE_immune_00001, prop_DE_t1d_00001]
categories = ['All Genes', 'Immune Genes', 'T1D Genes']

# Create subplots
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Plot for 0.001
ax[0].bar(categories, proportions_001, color='skyblue')
ax[0].set_title('Sex Effect - Proportions of DE Genes at 0.001 (Cutoff 0.2)')
ax[0].set_ylabel('Proportion of DE Genes')
ax[0].set_ylim(0, 0.025)

# Plot for 0.00001
ax[1].bar(categories, proportions_00001, color='lightgreen')
ax[1].set_title('Sex Effect - Proportions of DE Genes at 0.00001 (Cutoff 0.2)')
ax[1].set_ylabel('Proportion of DE Genes')
ax[1].set_ylim(0, 0.025)

plt.tight_layout()

#Save fig
plt.savefig(inDir+"prop_genes_DE_sex_effect.png")
plt.show()
plt.close()

