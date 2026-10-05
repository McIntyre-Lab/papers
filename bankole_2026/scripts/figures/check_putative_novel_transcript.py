#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 12:12:17 2026

@author: mgaran
"""
import pandas as pd

PROJ = "/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"

anno = pd.read_csv(f"{PROJ}/zenodo/fiveSpecies_dmel6_full_annotation.csv")
data = pd.read_csv(f"{PROJ}/zenodo/datafiles/datafile_jxnHash_dmel6.csv")
netanya = pd.read_csv(f"{PROJ}/sqanti_reads_analysis/erp_FBgn0034928_+11.csv")
netanya.columns.tolist()
['dmel_geneID',
 'dsim_geneID',
 'dyak_geneID',
 'yak_merge',
 'dsan_geneID',
 'san_merge',
 'dser_geneID',
 'ser_merge',
 'flag_novel_transcript_in_mel',
 'novel_jxnHash_mel',
 'flag_novel_transcript_in_sim',
 'novel_jxnHash_sim',
 'flag_novel_transcript_in_yak',
 'novel_jxnHash_yak',
 'flag_novel_transcript_in_san',
 'novel_jxnHash_san',
 'flag_novel_transcript_in_ser',
 'novel_jxnHash_ser',
 'flag_geneID_not_in_ortholog_list',
 'flag_any_species_novel',
 'flag_novel_in_gt1_species',
 'ERP_mel',
 'ERP_sim',
 'ERP_yak',
 'ERP_san',
 'ERP_ser']

# get putative novel jxnHash from netanya's squanti reads analysis output
jxnHash = netanya['novel_jxnHash_mel'].iloc[0] # only 1 row (minus header) so get first value in column
print(jxnHash)
# c89208844365baf911a4dfa42da435225a019490ec8f2ecec06129aab626daed

# get geneID of the jxnHash in novel jxnHash
gene1 = data[data['dmel6_jxnHash'] == jxnHash]['geneID'].iloc[0]
print(gene1)
# FBgn0267588

# check geneID in netanya's file
gene2 = netanya['dmel_geneID'].iloc[0] # only 1 row (minus header) so get first value in column
print(gene2)
# FBgn0034928

#### This gene is different then the gene assignment from netanyas sqanti reads analysis for the same jxnHash

# check geneset?
geneset1 = anno[anno['geneID'] == gene1]['genesetid'].iloc[0]
geneset2 = anno[anno['geneID'] == gene2]['genesetid'].iloc[0]
print(geneset1)
print(geneset2)

#### both geneID have geneset 23

# how many genes in geneset 23?

print(anno[anno['genesetid'] == geneset1]['geneID'].nunique())
# 4936

# get the annotated jxnHash in both genes: # FBgn0267588 and FBgn0034928
jxnHash_list = anno[anno['geneID'].isin([gene1, gene2])]["dmel6_jxnHash"].tolist()
print(jxnHash_list)

# get the components of the annotated jxnHash, convert to integers, and store in a list for plotting across species
component_list = anno[anno['dmel6_jxnHash'].isin(jxnHash_list)]['component_id'].astype(int).unique().tolist()
print(component_list)
[7295, 36119]

print(anno[anno['dmel6_jxnHash'] == 'f0e5c2cf4aa63da8f984da02b52dd71f9c4a834fff18a483da32f959ced217f9']['component_id'].iloc[0])