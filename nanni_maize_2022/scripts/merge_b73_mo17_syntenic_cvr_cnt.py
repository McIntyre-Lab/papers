#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Merge B73 and Mo17 coverage count summary values on syntenic list."
    )

    # Input data
    parser.add_argument(
        "-b",
        "--b73",
        dest="inB",
        required=True,
        help=(
            "Input CSV of B73 mapped gene-level coverage coun file in wide "
            "format and with summary values."
        )
    )
    parser.add_argument(
        "-m",
        "--mo17",
        dest="inM",
        required=True,
        help=(
            "Input CSV of Mo17 mapped gene-level coverage coun file in wide "
            "format and with summary values."
        )
    )
    parser.add_argument(
        "-s",
        "--synteny",
        dest="inSyn",
        required=True,
        help=(
            "B73-Mo17 CAU Synteny list from Nature Genetics (no header)."
        )
    )
    parser.add_argument(
        "-g",
        "--gene",
        dest="inGene",
        required=True,
        help=(
            "CSV file of genes of interest (column name gene_id)."
        )
    )
    parser.add_argument(
        "-t",
        "--treatment",
        dest="inTrt",
        required=True,
        choices=["Amb", "Ele"],
        help=(
            "Selected treatment to plot: Amb or Ele for ambient or elevated ozone."
        )
    )
    parser.add_argument(
        "-v",
        "--variable",
        dest="inV",
        required=True,
        default = "mean_tpm",
        help=(
            "Variable prefix for value to plot (Default: mean_tpm)."
        )
    )
    # Output data
    parser.add_argument(
        "-1",
        "--1to1",
        dest="out1to1",
        required=False,
        action="store_true",
        help=(
            "Output only 1-to-1 Mo17-B73 pairs (Default: output all pairs "
            "including some that may be paired with another gene)."
        )
    )
    parser.add_argument(
        "-p",
        "--prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix."
    )
    parser.add_argument(
        "-c",
        "--counts",
        dest="outCounts",
        required=False,
        help="Output CSV for counts of genes within 1%, 5%, or greater in eithere B73 or Mo17 mapping."
    )

    args = parser.parse_args()
    return args

def main():
    # Get input wide coverage count files
    inB = pd.read_csv(args.inB, low_memory=False)
    inM = pd.read_csv(args.inM, low_memory=False)
    treatment = args.inTrt
    valCol = args.inV

    # Get genes of interest
    geneDF = pd.read_csv(args.inGene, low_memory=False)[["gene_id"]]

    # Synteny list from Nature Genetics
    # Single column of gene ids, using gene as the column name to match AMM in B73v4_Mo17CAU_synteny_NatureGenetics.sas
#    synListNG = pd.read_csv(args.inSyn, names=["gene"])
    # ~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4/Mo17_Table2_gene_list/1.B73_Syntenic_genes/19.B73-sytenic_gene.txt
    # Add synteny flag
#    synListNG["flag_Mo17_B73_synteny_NatureGenetics"] = 1
    synDF = pd.read_csv(args.inSyn, sep="\t")

    # Combine synfnd and and synmap gene matches
    conditions = [
        synDF["match_synfind_synmap"] == "only synfind",
        synDF["match_synfind_synmap"] == "only synmap",
        synDF["match_synfind_synmap"] == "match"
    ]
    choices = [
        synDF["Mo17_CAU_Synfind_default"],
        synDF["Mo17_CAU_Synmap_megablast_0_001_"],
        synDF["Mo17_CAU_Synmap_megablast_0_001_"]
    ]
    synDF["Mo17_gene_id"] = np.select(conditions, choices, np.nan)

    # Count pairs that are 1-to-1 and drop those that are not if requested
    synDF["num_B73_to_Mo17"] = synDF.groupby("Mo17_gene_id")["B73v4_gene_ID"].transform("nunique")
    synDF["num_Mo17_to_B73"] = synDF.groupby("B73v4_gene_ID")["Mo17_gene_id"].transform("nunique")
    print("{} total B73-Mo17 syntenic pairs\n{} pairs are one-to-one".format(
        len(synDF[synDF["num_B73_to_Mo17"] + synDF["num_Mo17_to_B73"] > 0]),
        len(synDF[synDF["num_B73_to_Mo17"] + synDF["num_Mo17_to_B73"] == 2])
    ))
    if args.out1to1:
        print("!!! Selecting only one-to-one pairs")
        synDF = synDF[synDF["num_B73_to_Mo17"] + synDF["num_Mo17_to_B73"] == 2]
    else:
        print("!!! Plotting ALL pairs (including pairs not one-to-one)")

    # Merge genes of interest with synteny list
    geneSynMerge = pd.merge(
        synDF,
        geneDF,
        how = "outer",
        left_on = "B73v4_gene_ID",
        right_on = "gene_id",
        validate = "1:1",
        indicator = "merge_check"
    )
    print("{} of {} genes of interest in syntenic list".format(
        len(geneSynMerge[geneSynMerge["merge_check"]=="both"]),
        len(geneDF)
    ))
    #geneSynMerge["merge_check"].value_counts()
    # Select only syntenic pairs that are in gene list of interest
    geneSyn = geneSynMerge[geneSynMerge["merge_check"]=="both"].drop(columns=["merge_check"])

    # Merge B73 coverage counts to geneSyn
    geneSynBmerge = pd.merge(
        inB,
        geneSyn,
        how = "outer",
        left_on = "primary_FBgn",
        right_on = "gene_id",
        validate = "1:1",
        indicator = "merge_check"
    )
    print("{} of {} genes of interest in syntenic list are in B73 coverage counts".format(
        len(geneSynBmerge[geneSynBmerge["merge_check"]=="both"]),
        len(geneSynMerge[geneSynMerge["merge_check"]=="both"])
    ))
    #geneSynBmerge["merge_check"].value_counts()
    # Select only syntenic pairs that are in gene list of interest
    geneSynB = geneSynBmerge[geneSynBmerge["merge_check"]=="both"].drop(columns=["merge_check"])

    # Merge Mo17 coverage counts to geneSynB
    if args.out1to1:
        geneSynBMmerge = pd.merge(
            inM,
            geneSynB,
            how = "outer",
            left_on = "gene_id",
            right_on = "Mo17_gene_id",
            suffixes = ["_mapMo17", "_mapB73"],
            validate = "1:1",
            indicator = "merge_check"
        )
    else:
        geneSynBMmerge = pd.merge(
            inM,
            geneSynB,
            how = "outer",
            left_on = "gene_id",
            right_on = "Mo17_gene_id",
            suffixes = ["_mapMo17", "_mapB73"],
            validate = "1:m",
            indicator = "merge_check"
        )
    print("{} of {} genes of interest in syntenic list and B73 coverage counts are in Mo17 coverage counts".format(
        len(geneSynBMmerge[geneSynBMmerge["merge_check"]=="both"]),
        len(geneSynBmerge[geneSynBmerge["merge_check"]=="both"])
    ))
    #geneSynBMmerge["merge_check"].value_counts()
    # Select only syntenic pairs that are in gene list of interest
    geneSynBM = geneSynBMmerge[geneSynBMmerge["merge_check"]=="both"].drop(columns=["merge_check"])

    # Get ratio of B73 mapped over Mo17 mapped
    geneCols = [c for c in geneSynBM.columns if "gene_id_" in c]
    summaryCols = [c for c in geneSynBM.columns if valCol in c and treatment in c]
    B73Col = [c for c in summaryCols if "_mapB73" in c][0]
    Mo17Col = [c for c in summaryCols if "_mapMo17" in c][0]
    geneSynBM["ratio_" + treatment + "_" + valCol + "_B73_Mo17"] = geneSynBM[B73Col] / geneSynBM[Mo17Col]

    # Ouput mean columns for selected treatment
    geneSynBM[geneCols + summaryCols + ["ratio_" + treatment + "_" + valCol + "_B73_Mo17"]].to_csv(
        "{}_{}.csv".format(args.outPrefix, valCol), index=False)

    # Print descriptive values of gene summary values
    print("Descriptive values of gene values and ratio:\n{}".format(
        geneSynBM[summaryCols + ["ratio_" + treatment + "_" + valCol + "_B73_Mo17"]].describe().to_string()
    ))
    print("\nDistribution of ratios (B73 mapped / Mo17 mapped) where both are greater than 0:\n{}".format(
        geneSynBM[
            (geneSynBM[B73Col] > 0)
            & (geneSynBM[Mo17Col] > 0)]["ratio_" + treatment + "_" + valCol + "_B73_Mo17"].describe().to_string()
    ))

    if args.outCounts is not None:
        outCounts = pd.Series(index=[
            "num_gene_"+valCol+"_gt0_mapB73",
            "num_gene_"+valCol+"_gt0_mapMo17",
            "num_gene_"+valCol+"_gt0_both",
            "num_gene_"+valCol+"_1perc_both",
            "num_gene_"+valCol+"_5perc_both",
            "num_gene_"+valCol+"_mapB73_greater",
            "num_gene_"+valCol+"_mapMo17_greater"],
            dtype=int
        )
        outCounts["num_gene_"+valCol+"_gt0_mapB73"] = geneSynBM[geneSynBM[B73Col]>0]["B73v4_gene_ID"].nunique()
        outCounts["num_gene_"+valCol+"_gt0_mapMo17"] = geneSynBM[geneSynBM[Mo17Col]>0]["Mo17_gene_id"].nunique()
        outCounts["num_gene_"+valCol+"_gt0_both"] = len(geneSynBM[
                (geneSynBM[B73Col]>0)
                & (geneSynBM[Mo17Col]>0)
            ])
        geneSynBM["prop_diff"] = abs(geneSynBM[B73Col] - geneSynBM[Mo17Col]) / ((geneSynBM[B73Col] + geneSynBM[Mo17Col])/2)
        outCounts["num_gene_"+valCol+"_1perc_both"] = len(geneSynBM[geneSynBM["prop_diff"]<=0.01])
        outCounts["num_gene_"+valCol+"_5perc_both"] = len(geneSynBM[geneSynBM["prop_diff"]<=0.05])
        outCounts["num_gene_"+valCol+"_mapB73_greater"] = len(
                geneSynBM[
                    (geneSynBM["prop_diff"]>0.05)
                    & (geneSynBM[B73Col] > geneSynBM[Mo17Col])
        ])
        outCounts["num_gene_"+valCol+"_mapMo17_greater"] = len(
                geneSynBM[
                    (geneSynBM["prop_diff"]>0.05)
                    & (geneSynBM[B73Col] < geneSynBM[Mo17Col])
        ])
        pd.DataFrame(outCounts).T.to_csv(args.outCounts, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
