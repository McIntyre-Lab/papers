#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:23:08 2024

@author: nkeil
"""
import os
import pandas as pd

inDir="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType"
annotFile="allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv"
annot_DF=pd.read_csv(os.path.join(inDir, annotFile), sep=",", header = 0, low_memory=False)

pdiffFile_CD4="quantify_t1d_pacbio_transcripts/ES_SV_CD4_pdiff_4_ma.csv"
pdiff_DF_CD4=pd.read_csv(os.path.join(inDir, pdiffFile_CD4), sep=",", header = 0, low_memory=False)

pdiffFile_CD8="quantify_t1d_pacbio_transcripts/ES_SV_CD8_pdiff_4_ma.csv"
pdiff_DF_CD8=pd.read_csv(os.path.join(inDir, pdiffFile_CD8), sep=",", header = 0, low_memory=False)

## Stack pdiff DFs
pdiff_DF=pd.concat([pdiff_DF_CD4, pdiff_DF_CD8] , axis=0, ignore_index=True)

## Make geneID column
pdiff_DF['gene_id'] = pdiff_DF['featureID'].str.split(':').str[0]

##Merge in annotations from annotation file

flag_ir_DF= annot_DF[['ef_id', 'ef_ir_flag']]

pdiff_annot_DF=pdiff_DF.merge(flag_ir_DF, left_on='featureID', right_on='ef_id', how='left')

#Make datframe with counts of exons introns and the sum of both

# Keep one row per unique feature per gene per cell type
feature_count_input_DF = pdiff_annot_DF.drop_duplicates(
    subset=['gene_id', 'cellType', 'featureID', 'ef_ir_flag']
)

count_DF = feature_count_input_DF.groupby(['gene_id','cellType']).agg({
    'featureID': 'size',  # Counts unique features within each gene/cell type
    'ef_ir_flag': [
        lambda x: (x == 0).sum(), 
        lambda x: (x == 1).sum()   
    ]
}).reset_index()



count_DF.columns = ['gene_id', 'cellType', 'count_total_frags', 'count_exon_frags', 'count_intron_frags']

#Import t1d and immune gene lists
t1d_gene="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists/T1D_candidate_genes_robertson_2021_ensembl.txt"
immune_gene="/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/candidate_gene_lists/immunobase_hg19_ensembl_entrez.csv"

t1d_geneDF=pd.read_csv(t1d_gene, sep=',', header=0, low_memory=False)
immune_geneDF=pd.read_csv(immune_gene, sep=',', header=0, low_memory=False)

##Load file with pvals for all features
pval_allFile="quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary/meta_analysis_pvals_moderator_heterogenity.csv"
pval_allDF=pd.read_csv(os.path.join(inDir, pval_allFile), sep=",", header = 0, low_memory=False)

#Rename columns
pval_allDF.columns = ['gene_id', 'cellType', 'pval_het_overall', 'pval_mod_overall', 'flag_sig_het_overall', 'flag_sig_mod_overall']

##Load file with pvals for metanalysis with exons only
pval_exonFile="quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_exonsOnly/meta_analysis_pvals_moderator_heterogenity_exonsOnly.csv"
pval_exonDF=pd.read_csv(os.path.join(inDir, pval_exonFile), sep=",", header = 0, low_memory=False)

#Rename columns

pval_exonDF.columns = ['gene_id', 'cellType', 'pval_het_exon', 'pval_mod_exon', 'flag_sig_het_exon', 'flag_sig_mod_exon']

##Load file with pvals for metanalysi with exons only
pval_intronFile="quantify_t1d_pacbio_transcripts/meta_analysis_DE_pdiffs_summary_intronsOnly/meta_analysis_pvals_moderator_heterogenity_intronsOnly.csv"
pval_intronDF=pd.read_csv(os.path.join(inDir, pval_intronFile), sep=",", header = 0, low_memory=False)

#Rename columns
pval_intronDF.columns = ['gene_id', 'cellType', 'pval_het_intron', 'pval_mod_intron', 'flag_sig_het_intron', 'flag_sig_mod_intron']

#Combine all 3 dataframes 
merged_DF = pval_allDF.merge(pval_exonDF, on=['gene_id','cellType'], how='left').merge(pval_intronDF, on=['gene_id','cellType'], how='left')

##Merge in counts
merged_DF_w_count=merged_DF.merge(count_DF, on=['gene_id','cellType'], how='left')

#Make flag gene columns
merged_DF_w_count['flag_t1d_gene']=merged_DF_w_count['gene_id'].isin(t1d_geneDF['ensembl_geneID']).astype(int)
merged_DF_w_count['flag_immune_gene']=merged_DF_w_count['gene_id'].isin(immune_geneDF['Ensembl_ID']).astype(int)

CD4_DF= merged_DF_w_count[merged_DF_w_count['cellType'] == 'CD4']
CD8_DF= merged_DF_w_count[merged_DF_w_count['cellType'] == 'CD8']

CD4_DF.to_csv(os.path.join(inDir, "quantify_t1d_pacbio_transcripts/metaanalysis_summary_CD4.csv"), index=False)
CD8_DF.to_csv(os.path.join(inDir, "quantify_t1d_pacbio_transcripts/metaanalysis_summary_CD8.csv"), index=False)