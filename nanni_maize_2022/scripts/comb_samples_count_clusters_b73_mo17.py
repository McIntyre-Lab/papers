#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Combine samples and compare clustering mapping in B73-Mo17 CAU syntenic genes."
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input-direectory",
        dest="inDir",
        required=True,
        help="Input directory for files of clusters compared to B73-Mo17 CAU synteny for each sample."
    )

    args = parser.parse_args()
    return args

def list_files(directory, extension, allow_ext=False):
    if allow_ext == True:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension) or "."+extension+"." in f)
    else:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension))

def main():
    # Get input directory
    #DIR = "/nfshome/adalena.nanni/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/maize_ainsworth/compare_b73_2_mo17/compare_clusters_in_synteny"
    DIR = args.inDir

    # Loop over *_compare.csv files in input directory
    isFirst = True
    for file in list_files(DIR, "csv"):
        # Read in file and add ext based on file name
        name = "_".join(file.split("_")[0:3])
        inDF = pd.read_csv(DIR + "/" + file, low_memory=False)
        inDF["sample"] = name
        if isFirst:
            concatDF = inDF.copy()
            isFirst = False
            continue
        else:
            concatDF = pd.concat([concatDF, inDF], ignore_index=True)

    # Get counts on combined
    outFile = open(DIR + "/combined_counts.txt", 'w')
    outFile.write("Input Files:{}".format(
            list(list_files(DIR, "csv"))
        )
    )

    # Count total clusters and cluster in 12604
    outFile.write("{} total clusters\n"
                  "{} clusters in {} B73 annotated genes (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} Mo17 CAU annotated genes (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} B73 and Mo17 CAU annotated gene pairs (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} of the 12604 genes of interest\n\n".format(
            concatDF["id"].nunique(),
            concatDF[concatDF["flag_annotated_gene_b73"] == 1]["id"].nunique(),
            concatDF[concatDF["flag_annotated_gene_b73"] == 1]["associated_gene_b73"].nunique(),
            concatDF[concatDF["flag_annotated_gene_mo17CAU"] == 1]["id"].nunique(),
            concatDF[concatDF["flag_annotated_gene_mo17CAU"] == 1]["associated_gene_mo17CAU"].nunique(),
            concatDF[concatDF["flag_annotated_gene_b73"] + concatDF["flag_annotated_gene_mo17CAU"] == 2]["id"].nunique(),
            len(concatDF[concatDF["flag_annotated_gene_b73"] + concatDF["flag_annotated_gene_mo17CAU"] == 2][["associated_gene_b73", "associated_gene_mo17CAU"]].drop_duplicates()),
            concatDF[concatDF["merge_check2"]=="both"]["id"].nunique(),
            concatDF[concatDF["merge_check2"]=="both"]["gene_id"].nunique()
        )
    )

    # Get pbid assignment counts for all clusters and cluster in 12604
    outFile.write(
        "PBid assignment categories for all clusters:\n{}\n\n"
        "PBid assignment categories for clusters in the {} of 12604 genes of "
        "interest:\n{}\n\n".format(
            concatDF["pbid_assignment"].value_counts().to_string(),
            concatDF[concatDF["merge_check2"]=="both"]["gene_id"].nunique(),
            concatDF[concatDF["merge_check2"]=="both"]["pbid_assignment"].value_counts().to_string()
        )
    )

    # Count number of clusters and genes matching between B73 and M017 mapping
    outFile.write(
        "{} clusters are in {} matching B73-Mo17CAU syntenic genes\n"
        "{} clusters are in {} matching B73-Mo17CAU syntenic genes in the {} "
        "of 12604 genes of interest\n\n".format(
            int(concatDF[~concatDF["id"].isna()]["flag_Mo17_B73_synteny_synfind_synmap"].sum()),
            concatDF[
                    (~concatDF["id"].isna())
                    & (concatDF["flag_Mo17_B73_synteny_synfind_synmap"]==1)
                ]["B73v4_gene_ID"].nunique(),
            int(concatDF[
                    (~concatDF["id"].isna())
                    & (concatDF["merge_check2"]=="both")
                ]["flag_Mo17_B73_synteny_synfind_synmap"].sum()),
            concatDF[
                    (~concatDF["id"].isna())
                    & (concatDF["merge_check2"]=="both")
                    & (concatDF["flag_Mo17_B73_synteny_synfind_synmap"]==1)
                ]["B73v4_gene_ID"].nunique(),
            concatDF[concatDF["merge_check2"]=="both"]["gene_id"].nunique()
        )
    )

    outFile.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


