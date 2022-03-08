#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Subset file for B73-Mo17 1-to-1 syntenic genes."
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input CSV file."
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
        "--genotype",
        dest="inGeno",
        required=True,
        choices=["B73", "Mo17"],
        help=(
            "Genotype to merge B73-Mo17 syntenic list (B73 or Mo17)."
        )
    )
    parser.add_argument(
        "-c",
        "--column",
        dest="inCol",
        required=True,
        help=(
            "Column name in input file to merge B73-Mo17 syntenic list on."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output CSV file."
    )

    args = parser.parse_args()
    return args

def main():
    # Get 1-to-1 pairs of B73-Mo17 syntenic list
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
    synDF = synDF[synDF["num_B73_to_Mo17"] + synDF["num_Mo17_to_B73"] == 2]

    # Get input file
    inDF = pd.read_csv(args.inFile, low_memory=False)

    # Merge input file with 1-to-1 B73-Mo17 syntenic pairs
    if args.inGeno == "B73":
        geneSynMerge = pd.merge(
            synDF,
            inDF,
            how = "outer",
            left_on = "B73v4_gene_ID",
            right_on = args.inCol,
            validate = "1:1",
            indicator = "merge_check"
        )
    elif args.inGeno == "Mo17":
        geneSynMerge = pd.merge(
            synDF,
            inDF,
            how = "outer",
            left_on = "Mo17_gene_id",
            right_on = args.inCol,
            validate = "1:1",
            indicator = "merge_check"
        )
    # Select only genes in synteny list
    geneSyn = geneSynMerge[geneSynMerge["merge_check"]=="both"].drop(columns=["merge_check"])
    print("{} total genes in input file\n{} genes are in one-to-one syntenic".format(
        len(inDF),
        len(geneSyn)
    ))

    # Drop columns from synteny file and output to file
    geneSyn.drop(columns=synDF.columns).to_csv(args.outFile, index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
