#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Replace gene value with an associated transcript value in "
            "tappAS input expression matrix file."
        )
    )

    # Input data
    parser.add_argument(
        "-e",
        "--expression",
        dest="inExp",
        required=True,
        help=(
            "Input expression file to replace gene values with in first "
            "column and remove first column header (requried by tappAS)."
        )
    )
    parser.add_argument(
        "-g2t",
        "--gene2transcript",
        dest="inG2T",
        required=False,
        help=(
            "Input CSV file of associated gene (gene_id) and single transcript "
            "value (first_transcript_id)."
        )
    )
    parser.add_argument(
        "-g",
        "--gtf",
        dest="inGTF",
        required=False,
        help=(
            "Input GTF file to get associated gene and transcript values from "
            "(first transcript value after sorting by gene and transcript names)."
        )
    )
    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output file name for resulting expression TSV."
    )

    args = parser.parse_args()
    return args

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

def main():
    # Use input gtf to gene_id to first_transcript_id associations
    if args.inG2T is not None:
        gene2transcript = pd.read_csv(args.inG2T)
    elif args.inGTF is not None:
        gene2transcript = get_gene_from_gtf(args.inGTF)
    else:
        exit("!!!ERROR: Must provide a gene-to-transcript file or GTF file.")

    # Get input expression matrix
    expDF = pd.read_csv(args.inExp, sep="\t")

    # Merge by gene_id and PRIMARY_FBGN
    gene2xcprtExp = pd.merge(
            gene2transcript,
            expDF,
            how="outer",
            left_on="gene_id",
            right_on="PRIMARY_FBGN",
            indicator="merge_check",
            validate="1:1"
    )
    if len(expDF) != gene2xcprtExp["merge_check"].value_counts()["both"]:
        print(
            "!!!WARNING: Genes in expression file are missing in GTF..."
            "{} genes missing from GTF".format(
                    gene2xcprtExp["merge_check"].value_counts()["right_only"]
            )
        )
    gene2xcprtExp = gene2xcprtExp[gene2xcprtExp["merge_check"]=="both"]

    # Output expression file with gene value replaced with first transcript
    #   and the header of the first column empty (required by tappAS)
    expCols = [c for c in expDF.columns if c != "PRIMARY_FBGN"]
    outExp = gene2xcprtExp[["first_transcript_id"] + expCols].rename(
        columns={
            "first_transcript_id": ""
        }
    )
    outExp.to_csv(args.outFile, sep="\t", index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
