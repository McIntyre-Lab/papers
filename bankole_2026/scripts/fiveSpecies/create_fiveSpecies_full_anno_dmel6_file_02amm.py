#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:15:06 2024

@author: ammorse
"""

import argparse
import pandas as pd

# Initialize argument parser
parser = argparse.ArgumentParser(description="Create fiveSpecies dmel6 full annotation.")

# Define command-line arguments
parser.add_argument('--jxnhash_fbtr_file', required=True, 
                    help="Path to the jxnHash to FBTR link file")
parser.add_argument('--genome', required=True, 
                    help="Name of the genome")
parser.add_argument('--erp_file', required=True, 
                    help="Path to the ERP file")
parser.add_argument('--esp_file', required=True, 
                    help="Path to the ESP file")
parser.add_argument('--flag_file', required=True, 
                    help="Path to the flag file")
parser.add_argument('--output_file', required=True, 
                    help="Path to save the output file")
parser.add_argument('--mismatch', required=True, 
                    help="Path to save the geneID mismatch file")

# Parse arguments
args = parser.parse_args()

# Access the arguments using args.<argument_name>
jxnhash_fbtr_file = args.jxnhash_fbtr_file
genome = args.genome
erp_file = args.erp_file
esp_file = args.esp_file
flag_file = args.flag_file
output_file = args.output_file
output_mismatch = args.mismatch

#genome="dmel6"
#x="dmel650"
#y="dmel6"
#jxnhash_fbtr_file = f"/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/{x}_2_{y}_ujc_xscript_link.csv" 
#erp_file =    f"/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/fiveSpecies_2_{y}_ujc_er_vs_fiveSpecies_2_{y}_ujc_infoERP.csv" 
#esp_file =   f"/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/fiveSpecies_2_{y}_ujc_es_vs_fiveSpecies_2_{y}_ujc_infoESP.csv" 
#flag_file =   f"/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/flag_fiveSpecies_2_{y}_ujc_w_IR_flag.csv" 


print(f"RUNNING {genome}")

# Import the link between the jxnHash and FBTR
link_jxnhash2fbtr = pd.read_csv(jxnhash_fbtr_file)

# Rename and drop columns
newNameHash = f"{genome}_jxnHash"
newNameGene = f"{genome}_geneID_anno"
link_jxnhash2fbtr.rename(columns={
    'jxnHash': newNameHash,
    'geneID': newNameGene
}, inplace=True)
#link_jxnhash2fbtr.drop(columns=['jxnString'], inplace=True)

# Sort and cat
link_jxnhash2fbtr_sorted = link_jxnhash2fbtr.sort_values(by=[newNameHash, "transcriptID"])
link_jxnhash2fbtr_sorted = link_jxnhash2fbtr_sorted.groupby([
    f"{genome}_jxnHash",f"{genome}_geneID_anno" ])["transcriptID"].agg(list).reset_index()
link_jxnhash2fbtr_sorted[f"cat_{genome}_transcriptID"] = link_jxnhash2fbtr_sorted["transcriptID"].apply(lambda x: '|'.join(x))

# drop transcriptID col
link_jxnhash2fbtr_sorted.drop(columns=['transcriptID'], inplace = True)

# check file is uniq on genome_jxnHash
is_unique = link_jxnhash2fbtr_sorted[f"{genome}_jxnHash"].is_unique
print("is unique?", is_unique)

# Find any non-unique values for jxnHash - geneID
non_unique_values = link_jxnhash2fbtr_sorted[link_jxnhash2fbtr_sorted[f"{genome}_jxnHash"].isin(link_jxnhash2fbtr_sorted[f"{genome}_jxnHash"].value_counts()[link_jxnhash2fbtr_sorted[f"{genome}_jxnHash"].value_counts() > 1].index)]
print("num of nonunique is: ", len(non_unique_values))

# Count ujc 
print(f"num uniq ujc in {genome} xscript link file: ",len(link_jxnhash2fbtr_sorted))    ## 30508

# Import the five species information
link_jxnhash2ERP = pd.read_csv(erp_file)
link_jxnhash2ESP = pd.read_csv(esp_file)
flag_fivespecies = pd.read_csv(flag_file)

# Sort the ERP and ESP dataframes by jxnHash
link_jxnhash2ERP.sort_values(by='jxnHash', inplace=True)
link_jxnhash2ESP.sort_values(by='jxnHash', inplace=True)

# Merge ERP and ESP links on jxnHash
jxn_annot = pd.merge(
    link_jxnhash2ERP, 
    link_jxnhash2ESP, 
    on='jxnHash', 
    how='outer',
    suffixes=('_ERP', '_ESP'),
    indicator='merge1'
)
jxn_annot['merge1'].value_counts(dropna=False).sort_index()

# Rename column
jxn_annot.rename(columns={'jxnHash': newNameHash}, inplace=True)

# Sort flag data and FBTR data by dmel6_jxnHash
flag_fivespecies.sort_values(by=f"{genome}_jxnHash", inplace=True)
link_jxnhash2fbtr_sorted.sort_values(by=f"{genome}_jxnHash", inplace=True)

# Merge flag, annotated jxns, and FBTR data
merged_df = flag_fivespecies.merge(
    jxn_annot, on=f"{genome}_jxnHash", how='outer', indicator='merge2'
).merge(
    link_jxnhash2fbtr_sorted, on=f"{genome}_jxnHash", how='outer', indicator='merge3'
)

print(merged_df['merge1'].value_counts(dropna=False).sort_index())
print(merged_df['merge2'].value_counts(dropna=False).sort_index())
print(merged_df['merge3'].value_counts(dropna=False).sort_index())

# check that geneID == geneID_ERP == geneID_ESP == dmel6_geneID_anno
cols_2_compare = ['geneID', 'geneID_ERP', 'geneID_ESP', 'dmel6_geneID_anno']
merged_df['geneID_check'] = merged_df[cols_2_compare].apply(
    lambda row: row.dropna().nunique() <= 1, axis=1)

mismatch= merged_df[~merged_df['geneID_check']]
print("num mismatched geneIDs: ",len(mismatch))
# subset mismatch to only geneID columns for checking
mismatch_geneID=mismatch.loc[:, mismatch.columns.str.contains('geneID')]                               

# Drop merge check columns     
merged_df.drop(columns=['merge1', 'merge2', 'merge3'], inplace=True)

# Drop extra dmel geneID columns
merged_df.drop(columns=['geneID_ERP', 'geneID_ESP', f"{genome}_geneID_anno", 
                        'geneID_check'], inplace=True)

## count rows
print(f"number of rows in full {genome} anno:  {len(merged_df)}")      # 75986
print(f"number of rows in mismatch {genome}:  {len(mismatch_geneID)}") # 2
 
print(merged_df.columns)
## rename cols for LMM
merged_df = merged_df.rename(columns={
    "flagDataOnlyExon_ERP": "flagBonusExon_ERP",
    "flagDataOnlyExon_ESP": "flagBonusExon_ESP",
    "numDataOnlyExon_ERP": "numBonusExon_ERP",
    "numDataOnlyExon_ESP": "numBonusExon_ESP",
    "dataOnlyER_ID": "bonusER_ID",
    "dataOnlyES_ID": "bonusES_ID"
    })

#print(merged_df.columns)
# Save the results
merged_df.to_csv(output_file, index=False)
mismatch_geneID.to_csv(output_mismatch, index=False)

print("Processing complete. Output saved to:", output_file)



