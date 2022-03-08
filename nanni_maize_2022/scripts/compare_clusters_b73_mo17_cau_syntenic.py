#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Compare clustering mapping in B73-Mo17 CAU syntenic genes."
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input CSV file with gene_id columns of genes to subset counts for."
    )
    parser.add_argument(
        "-tB73",
        "--tofu-B73",
        dest="T_B73",
        required=True,
        help="Output file of tofu2 collapse (*.read_stat.txt) when mapped to B73 genome."
    )
    parser.add_argument(
        "-tMo17",
        "--tofu-Mo17",
        dest="T_Mo17",
        required=True,
        help="Output file of tofu2 collapse (*.read_stat.txt) when mapped to Mo17 CAU genome."
    )
    parser.add_argument(
        "-sB73",
        "--sqanti-B73",
        dest="S_B73",
        required=True,
        help="Output file of SQANTI QC (*_classification.txt) when mapped to B73 genome."
    )
    parser.add_argument(
        "-sMo17",
        "--sqanti-Mo17",
        dest="S_Mo17",
        required=True,
        help="Output file of SQANTI QC (*_classification.txt) when mapped to Mo17 CAU genome."
    )
    parser.add_argument(
        "--synteny",
        dest="Syn",
        required=True,
        help="TSV file from COGE B73v4-Mo17CAU synteny."
    )

    # Output data
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix for counts and output files."
    )

    args = parser.parse_args()
    return args

def main():
    # Get intput gene file
    geneDF = pd.read_csv(args.inFile)
#    geneDF = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/resubmission_2021/supplement/Supplementary_File_5.csv")

    # Synteny list from COGE
    synListC = pd.read_csv(args.Syn, sep="\t")
#    synListC = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv", sep="\t")
    # Add synteny flag
    synListC["flag_Mo17_B73_synteny_synfind_synmap"] = np.where(
            synListC["match_synfind_synmap"].isin(["match", "only synfind", "only synmatch"]),
            1,
            0
    )
    # Add Mo17 CAU gene id column (either synfind or synmap gene id)
    synListC["Mo17_CAU_gene_ID"] = np.where(
        synListC["match_synfind_synmap"].isin(["match", "only synfind"]),
        synListC["Mo17_CAU_Synfind_default"],
        np.where(
                synListC["match_synfind_synmap"] == "only synmap",
                synListC["Mo17_CAU_Synmap_megablast_0_001_"],
                np.nan
        )
    )

    # Tofu2 collapse read_stat files
    tB73 = pd.read_csv(args.T_B73, sep="\t")
    tMo17 = pd.read_csv(args.T_Mo17, sep="\t")
#    tB73 = pd.read_csv("/nfshome/adalena.nanni/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/maize_ainsworth/compare_b73_2_mo17/tofu2_b73/19_mo17_amb/19_mo17_amb.collapsed.read_stat.txt", sep = "\t")
#    tMo17 = pd.read_csv("/nfshome/adalena.nanni/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/maize_ainsworth/compare_b73_2_mo17/tofu2_mo17_cau/19_mo17_amb/19_mo17_amb.collapsed.read_stat.txt", sep = "\t")

    # SQANTI QC classification files
    sB73 = pd.read_csv(args.S_B73, sep = "\t")
    sMo17 = pd.read_csv(args.S_Mo17, sep = "\t")
#    sB73 = pd.read_csv("/nfshome/adalena.nanni/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/maize_ainsworth/compare_b73_2_mo17/sqanti_b73/mo17/amb/19_mo17_amb_classification.txt", sep = "\t")
#    sMo17 = pd.read_csv("/nfshome/adalena.nanni/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/maize_ainsworth/compare_b73_2_mo17/sqanti_mo17_cau/mo17/amb/19_mo17_amb_classification.txt", sep = "\t")

    # Add flag for annotation genes
    #   (structural_category is full-splice_match, incomplete-splice_match, novel_in_catalog, or novel_not_in_catalog)
    sB73["flag_annotated_gene"] = np.where(
        sB73["structural_category"].isin([
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ]),
        1,
        0
    )
    sMo17["flag_annotated_gene"] = np.where(
        sMo17["structural_category"].isin([
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ]),
        1,
        0
    )
    # For each of the 3 Mo17 samples, count clusters that do and do not map to matching Mo17-B73 syntenic genes using the following method:
    # 1) Use tofu2 read_stat file to associate raw long read "movie" ids to
    #       PBids for samples mapped to B73 and Mo17 CAU
    # 2) Merge in SQANTI QC classification file to associate PBids to reference
    #       gene_ids for samples mapped to B73 and Mo17 CAU
    tofuSqantiB73merge = pd.merge(
            tB73,
            sB73,
            how = "outer",
            left_on = "pbid",
            right_on = "isoform",
            validate = "m:1",
            indicator = "merge_check"
    )
    tofuSqantiMo17merge = pd.merge(
            tMo17,
            sMo17,
            how = "outer",
            left_on = "pbid",
            right_on = "isoform",
            validate = "m:1",
            indicator = "merge_check"
    )
#    tofuSqantiB73merge["merge_check"].value_counts()
#    both          240128
#    left_only      14675
#    right_only         0
#    tofuSqantiMo17merge["merge_check"].value_counts()
#    both          248673
#    left_only       6130
#    right_only         0

    # Count the number of movie ids per pbid
    tofuSqantiB73merge["num_id_per_pbid"] = tofuSqantiB73merge.groupby("pbid")["id"].transform("nunique")
    tofuSqantiMo17merge["num_id_per_pbid"] = tofuSqantiMo17merge.groupby("pbid")["id"].transform("nunique")

    # 3) Merge B73 and Mo17 CAU mapped variables by raw long read movie id
    #   Select only id (movie id), pdid, and associated_gene columns
    keepCols = ["id", "pbid", "associated_gene", "num_id_per_pbid", "flag_annotated_gene"]
    tofuSqantiB73Mo17Merge = pd.merge(
        tofuSqantiB73merge[keepCols],
        tofuSqantiMo17merge[keepCols],
        how = "outer",
        on = "id",
        suffixes = ["_b73", "_mo17CAU"],
        validate = "1:1",
        indicator = "merge_check"
    )
#    tofuSqantiB73Mo17Merge["merge_check"].value_counts()
#    both          254803
#    right_only         0
#    left_only          0

    # Add category for pbid assignment
    pbidCond = [
        (~tofuSqantiB73Mo17Merge["pbid_b73"].isna()) & (~tofuSqantiB73Mo17Merge["pbid_mo17CAU"].isna()),
        (~tofuSqantiB73Mo17Merge["pbid_b73"].isna()) & (tofuSqantiB73Mo17Merge["pbid_mo17CAU"].isna()),
        (tofuSqantiB73Mo17Merge["pbid_b73"].isna()) & (~tofuSqantiB73Mo17Merge["pbid_mo17CAU"].isna()),
        (tofuSqantiB73Mo17Merge["pbid_b73"].isna()) & (tofuSqantiB73Mo17Merge["pbid_mo17CAU"].isna())
    ]
    pbidAssign = ["both", "b73", "mo17CAU", "none"]
    tofuSqantiB73Mo17Merge["pbid_assignment"] = np.select(pbidCond, pbidAssign, "oops")
#    tofuSqantiB73Mo17Merge["pbid_assignment"].value_counts()
#    both       236911
#    mo17CAU     11762
#    b73          3217
#    none         2913

    # 4) Merge in B73_Mo17CAU file from COGE to associate B73 reference and
    #       Mo17 CAU reference gene_ids
    tofuSqantiB73Mo17Syn = pd.merge(
        tofuSqantiB73Mo17Merge.drop(columns=["merge_check"]),
        synListC,
        how = "outer",
        left_on = ["associated_gene_b73", "associated_gene_mo17CAU"],
        right_on = ["B73v4_gene_ID", "Mo17_CAU_gene_ID"],
        validate = "m:1",
        indicator = "merge_check"
    )
#    tofuSqantiB73Mo17Syn["merge_check"].value_counts()
#    both          171356
#    left_only      83447
#    right_only     42511

    # 5) Merge in 12604 genes on interest
    tofuSqantiB73Mo17Syn2 = pd.merge(
        tofuSqantiB73Mo17Syn,
        geneDF,
        how = "outer",
        left_on = "associated_gene_b73",
        right_on = "gene_id",
        validate = "m:1",
        indicator = "merge_check2"
    )
#    tofuSqantiB73Mo17Syn2["merge_check2"].value_counts()
#    both          216102
#    left_only      81212
#    right_only      6028

    # Output file for sample
    tofuSqantiB73Mo17Syn2.to_csv(args.outPrefix+".csv", index = False)

    # 6) Get counts
    outFile = open(args.outPrefix+"_counts.txt", 'w')
    outFile.write("Input Files:\n\t{}\n\t{}\n\t{}\n\t{}\n\t{}\n\t{}\n\n".format(
            args.inFile,
            args.Syn,
            args.T_B73,
            args.T_Mo17,
            args.S_B73,
            args.S_Mo17
        )
    )

    # Count total clusters and cluster in 12604
    outFile.write("{} total clusters\n"
                  "{} clusters in {} B73 annotated genes (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} Mo17 CAU annotated genes (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} B73 and Mo17 CAU annotated gene pairs (SQANTI QC FSM, ISM, NIC, or NNC)\n"
                  "{} clusters in {} of the 12604 genes of interest\n\n".format(
            tofuSqantiB73Mo17Syn2["id"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_b73"] == 1]["id"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_b73"] == 1]["associated_gene_b73"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_mo17CAU"] == 1]["id"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_mo17CAU"] == 1]["associated_gene_mo17CAU"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_b73"] + tofuSqantiB73Mo17Syn2["flag_annotated_gene_mo17CAU"] == 2]["id"].nunique(),
            len(tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["flag_annotated_gene_b73"] + tofuSqantiB73Mo17Syn2["flag_annotated_gene_mo17CAU"] == 2][["associated_gene_b73", "associated_gene_mo17CAU"]].drop_duplicates()),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["merge_check2"]=="both"]["id"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["merge_check2"]=="both"]["gene_id"].nunique()
        )
    )

    # Get pbid assignment counts for all clusters and cluster in 12604
    outFile.write(
        "PBid assignment categories for all clusters:\n{}\n\n"
        "PBid assignment categories for clusters in the {} of 12604 genes of "
        "interest:\n{}\n\n".format(
            tofuSqantiB73Mo17Syn2["pbid_assignment"].value_counts().to_string(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["merge_check2"]=="both"]["gene_id"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["merge_check2"]=="both"]["pbid_assignment"].value_counts().to_string()
        )
    )

    # Count number of clusters and genes matching between B73 and M017 mapping
    outFile.write(
        "{} clusters are in {} matching B73-Mo17CAU syntenic genes\n"
        "{} clusters are in {} matching B73-Mo17CAU syntenic genes in the {} "
        "of 12604 genes of interest\n\n".format(
            int(tofuSqantiB73Mo17Syn2[~tofuSqantiB73Mo17Syn2["id"].isna()]["flag_Mo17_B73_synteny_synfind_synmap"].sum()),
            tofuSqantiB73Mo17Syn2[
                    (~tofuSqantiB73Mo17Syn2["id"].isna())
                    & (tofuSqantiB73Mo17Syn2["flag_Mo17_B73_synteny_synfind_synmap"]==1)
                ]["B73v4_gene_ID"].nunique(),
            int(tofuSqantiB73Mo17Syn2[
                    (~tofuSqantiB73Mo17Syn2["id"].isna())
                    & (tofuSqantiB73Mo17Syn2["merge_check2"]=="both")
                ]["flag_Mo17_B73_synteny_synfind_synmap"].sum()),
            tofuSqantiB73Mo17Syn2[
                    (~tofuSqantiB73Mo17Syn2["id"].isna())
                    & (tofuSqantiB73Mo17Syn2["merge_check2"]=="both")
                    & (tofuSqantiB73Mo17Syn2["flag_Mo17_B73_synteny_synfind_synmap"]==1)
                ]["B73v4_gene_ID"].nunique(),
            tofuSqantiB73Mo17Syn2[tofuSqantiB73Mo17Syn2["merge_check2"]=="both"]["gene_id"].nunique()
        )
    )

    outFile.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
