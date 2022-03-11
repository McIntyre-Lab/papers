#!/usr/bin/env python3

# Get union of B73-Mo17 synteny lists and compare to 12604

import pandas as pd
import numpy as np

# Get input files
geneDF = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/resubmission_2021/supplement/Supplementary_File_5.csv")

# Synteny list from COGE
synListC = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv",
                       sep="\t")
# Add synteny flag
synListC["flag_Mo17_B73_synteny_synfind_synmap"] = np.where(
        synListC["match_synfind_synmap"].isin(["match", "only synfind", "only synmatch"]),
        1,
        0
)

# Synteny list from Nature Genetics
# Single column of gene ids, using gene as the column name to match AMM in B73v4_Mo17CAU_synteny_NatureGenetics.sas
synListNG = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt",
                       names=["gene"])
# Add synteny flag
synListNG["flag_Mo17_B73_synteny_NatureGenetics"] = 1

# Merge synteny lists
synMerge = pd.merge(
        synListC,
        synListNG,
        how = "outer",
        left_on = "B73v4_gene_ID",
        right_on = "gene",
        validate = "1:1",
        indicator = "merge_check"
)
synMerge["merge_check"].value_counts()
#     both          33681
#    left_only     15521
#    right_only        0

pd.crosstab(
        synMerge["flag_Mo17_B73_synteny_synfind_synmap"],
        synMerge["flag_Mo17_B73_synteny_NatureGenetics"],
        dropna = False
)
#    flag_Mo17_B73_synteny_NatureGenetics    1.0
#    flag_Mo17_B73_synteny_synfind_synmap       
#    0                                      5853
#    1                                     27828

# All syntenic genes from COGE is a subset of genes from the nature genetics list

# Using only nature genetics list moving forward

# Merge synteny list with 12604 and count
geneSynMerge = pd.merge(
        synListNG,
        geneDF,
        how = "outer",
        left_on = "gene",
        right_on = "gene_id",
        validate = "1:1",
        indicator = "merge_check"
)
geneSynMerge["merge_check"].value_counts()
#    left_only     21803
#    both          11878
#    right_only      726

# Drop the left_only, or genes only in the synteny list
geneSyn = geneSynMerge[geneSynMerge["merge_check"]!="left_only"].drop(columns=["merge_check"])

# Count genes in 12604 that are in synteny list (should be same as the "both")
geneSyn["flag_Mo17_B73_synteny_NatureGenetics"].sum()
# 11878
# (94% Mo17-B73 syntenic)
