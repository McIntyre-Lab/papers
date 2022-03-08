#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

def get_gene_from_gtf(infile):
    # Get input GTF file
    gtf = pd.read_csv(
            infile,
            names=[
                "chr",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute"
            ],
            comment="#",
            dtype=str,
            compression="infer",
            sep="\t",
            low_memory=False
    )
    # Select only exon features
    gtf = gtf[gtf["feature"]=="exon"]
    # Get attribute values
    gtf["attribute_values"] = gtf["attribute"].str.split(" ")
    # Check that gene_id and transcript_id attributes are present for all exons
    if not gtf["attribute_values"].apply(lambda x: "transcript_id" in x).all():
        print(
            "ERROR: transcript_id not contained within all exon "
            "attributes".format(infile)
        )
    if not gtf["attribute_values"].apply(lambda x: "gene_id" in x).all():
        print(
            "ERROR: gene_id not contained within all exon "
            "attributes".format(infile)
        )
    # Extract transcript_id and gene_id from attribute column
    gtf["transcript_id"] = gtf.apply(
            lambda x: x["attribute_values"][
                    x["attribute_values"].index("transcript_id")+1
                    ].split(";")[0].split("\"")[1],
            axis=1
    )
    gtf["gene_id"] = gtf.apply(
            lambda x: x["attribute_values"][
                    x["attribute_values"].index("gene_id")+1
                    ].split(";")[0].split("\"")[1],
            axis=1
    )
    # sort by gene and transcript id
    gtf = gtf.sort_values(["gene_id", "transcript_id"])
    # select first transcript_id per gene
    geneDF = gtf.groupby("gene_id")[
        [
            "chr",
            "transcript_id"
        ]].first().reset_index().rename(columns={
            "transcript_id": "first_transcript_id"
    })
    # Return unique list of genes with corresponding chromosome and first transcript
    return geneDF

# Get directories
inref = "/nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/useful_maize_info"
indata = "/nfshome/adalena.nanni/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/make_files_4_tappas_amm"

# Compare genes in maize B73 RefGen_v4.41 to v4.34
gtf41 = get_gene_from_gtf(inref+"/Zea_mays.B73_RefGen_v4.41.gtf.gz")
gtf34 = get_gene_from_gtf(inref+"/Zea_mays.AGPv4.34.gtf.gz")

# Output gene to first transcript CSV file for each annotation
gtf41[["gene_id", "first_transcript_id"]].to_csv(
        inref+"/Zea_mays.B73_RefGen_v4.41_gene_2_first_transcript.csv")
gtf34[["gene_id", "first_transcript_id"]].to_csv(
        inref+"/Zea_mays.AGPv4.34_gene_2_first_transcript.csv")

# Merge 41 and 34
gene_41_34 = pd.merge(
        gtf41,
        gtf34,
        how="outer",
        on="gene_id",
        suffixes=["_41", "_34"],
        indicator="merge_check",
        validate="1:1"
)
gene_41_34["merge_check"].value_counts()
#    both          44300
#    left_only      2130
#    right_only      174

gene_41_34["flag_in_41_only"] = np.where(
        gene_41_34["merge_check"]=="left_only",
        1,
        0
)
gene_41_34["flag_in_34_only"] = np.where(
        gene_41_34["merge_check"]=="right_only",
        1,
        0
)
gene_41_34["flag_in_both_41_34"] = np.where(
        gene_41_34["merge_check"]=="both",
        1,
        0
)
gene_41_34 = gene_41_34.drop(columns=["merge_check"])
gene_41_34[[c for c in gene_41_34.columns if "flag" in c]].sum()
#    flag_in_41_only        2130
#    flag_in_34_only         174
#    flag_in_both_41_34    44300

# Get list of genes going into tappas to merge with the 41 and 34 set
def list_files(directory, suffix, allow_ext=False, prefix=None):
    if prefix is None:
        return (f for f in os.listdir(directory) if f.endswith(suffix))
    else:
        return (f for f in os.listdir(directory) if (f.endswith(suffix)) and 
                (f.startswith(prefix)))
dataGeneDF = pd.DataFrame(columns=["gene_id"])
for file in list_files(directory=indata, suffix="_4_tappas_tpm.tsv", prefix="sbys"):
    # Get expression matrix and concat gene list
    tempDF = pd.read_csv(indata+"/"+file, sep="\t").rename(columns={
            "PRIMARY_FBGN": "gene_id"})
    dataGeneDF = pd.concat(
            [dataGeneDF, tempDF[["gene_id"]]],
            ignore_index=True,
            sort=True
    )
dataGeneDF = dataGeneDF.drop_duplicates()

# Merge tappas input genes with 41 and 34 set
dataGeneMerge4134 = pd.merge(
        gene_41_34,
        dataGeneDF,
        how="outer",
        on="gene_id",
        indicator="merge_check",
        validate="1:1"
)
dataGeneMerge4134["merge_check"].value_counts()
#    left_only     36751
#    both           9853
#    right_only        0

dataGeneMerge4134["flag_in_tappas"] = np.where(
        dataGeneMerge4134["merge_check"]=="both",
        1,
        0
)

print(
  dataGeneMerge4134.groupby(
          [c for c in dataGeneMerge4134.columns if "flag" in c]
  )["gene_id"].count().reset_index().to_string(index=False)
)
# flag_in_41_only  flag_in_34_only  flag_in_both_41_34  flag_in_tappas  gene_id
#               0                0                   1               0    34662
#               0                0                   1               1     9638
#               0                1                   0               0      174
#               1                0                   0               0     1915
#               1                0                   0               1      215

# 215 genes in the tappas files from 41 that are NOT in 34

