#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Plot transcript lengths from FASTA/FASTQ/CSV of lengths."
    )

    # Input data
    parser.add_argument(
        "-fq",
        "--fastq",
        dest="inFQ",
        action="append",
        required=False,
        help=(
            "Input FQ file(s) of interest, for each additional -fq argument "
            "the FQ lengths will be added to the same distribution plot (each "
            "will also get individual histograms)."
        )
    )
    parser.add_argument(
        "-fa",
        "--fasta",
        dest="inFA",
        action="append",
        required=False,
        help=(
            "Input FA file(s) of interest, for each additional -fa argument "
            "the FA lengths will be added to the same distribution plot (each "
            "will also get individual histograms)."
        )
    )
    parser.add_argument(
        "-l",
        "--length",
        dest="inLen",
        action="append",
        required=False,
        help=(
            "Input CSV file(s) with 'name' and 'length' columns, for each "
            "additional -l argument the lengths will be added to the same "
            "distribution plot (each will also get individual histograms)."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output-directory",
        dest="outDir",
        required=True,
        help="Output directory."
    )
    parser.add_argument(
        "-pdf",
        "--output-pdf",
        dest="outPDF",
        required=False,
        action="store_true",
        help="Optionally output the combined figure as botha PDF and PNG (Default: PNG only)."
    )
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=False,
        help=(
            "Output prefix for the combined distribution plot (not required - "
            "named all_seq_length_dist.png otherwise)."
        )
    )

    args = parser.parse_args()
    return args


def plot_hist(df, col, title=None):
    series = df[col]
    p = series.hist(bins=100, figsize=(10,4))
    plt.text(0.4, 0.95, "Mean = {}{}{}\nQ1 = {}\nMedian = {}\nQ3 = {}\nMax = {}".format(
            round(series.mean(), 3), u'\u00B1', round(series.std(), 3),
            round(series.quantile(.25), 3),
            round(series.median(), 3),
            round(series.quantile(.75), 3),
            df[col].max()), ha='left', va='top', transform=p.transAxes)
    if title is not None:
        plt.title(title)
    plt.tight_layout()

def plot_dist(df, col, label, hist=True):
    if sns.__version__ == "0.9.0":
        sns.distplot(
            df[col],
            hist=hist,
            label=label
        )
    else:
        print("WARNING: Script currently requires seaborn v0.9.0 to plot combined distribution...skipping.")

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def is_fastq(filename):
    with open(filename, "r") as handle:
        fastq = SeqIO.parse(handle, "fastq")
        try:
            return any(fastq)
        except Exception as e:
            print(e)
            return False

def main():
    # Get intput files
    fqList = args.inFQ
    faList = args.inFA
    lenList = args.inLen
    lengthDict = dict()

    # Process fastq files
    if fqList is not None:
        for fastq in fqList:
            # Check if properly formatted FA
            if not is_fastq(fastq):
                print("WARNING: File {} not properly formatted FASTQ...skipping.".format(fastq))
                continue
            # Read in FQ sequences and store lengths
            name = fastq.split("/")[-1].split(".")[0]
            lengthDF = pd.DataFrame(columns=["name", "length"])
            for record in SeqIO.parse(fastq, "fastq"):
                lengthDF = lengthDF.append(
                    pd.DataFrame(
                            [[record.id, len(record.seq)]], columns=["name", "length"]
                    ),
                    ignore_index=True
                )
            lengthDF.to_csv("{}/{}_seq_length.csv".format(args.outDir, name), index=False)
            lengthDict[name] = lengthDF
            # Plot sequence length distribution histogram for individual fastq
            plot_hist(lengthDF, "length", "Sequence Length Distribution")
            plt.savefig("{}/{}_seq_length_hist.png".format(args.outDir, name),dpi=600,format="png")
            plt.close()

    # Process fasta files
    if faList is not None:
        for fasta in faList:
            # Check if properly formatted FA
            if not is_fasta(fasta):
                print("WARNING: File {} not properly formatted FASTA...skipping.".format(fasta))
                continue
            # Read in FA sequences and store lengths
            name = fasta.split("/")[-1].split(".")[0]
            lengthDF = pd.DataFrame(columns=["name", "length"])
            for record in SeqIO.parse(fasta, "fasta"):
                lengthDF = lengthDF.append(
                    pd.DataFrame(
                            [[record.id, len(record.seq)]], columns=["name", "length"]
                    ),
                    ignore_index=True
                )
            lengthDF.to_csv("{}/{}_seq_length.csv".format(args.outDir, name), index=False)
            lengthDict[name] = lengthDF
            # Plot sequence length distribution histogram for individual fasta
            plot_hist(lengthDF, "length", "Sequence Length Distribution")
            plt.savefig("{}/{}_seq_length_hist.png".format(args.outDir, name),dpi=600,format="png")
            plt.close()

    # Process length list files
    if lenList is not None:
        for length in lenList:
            # Read in CSV and check that "length" column is present
            name = length.split("/")[-1].split(".")[0]
            lengthDF = pd.read_csv(length, low_memory=False)
            if "length" not in lengthDF.columns:
                print("WARNING: File {} does not contain 'length' column name...skipping.".format(length))
                continue
            else:
                lengthDict[name] = lengthDF
            # Plot sequence length distribution histogram for individual file
            plot_hist(lengthDF, "length", "Sequence Length Distribution")
            plt.savefig("{}/{}_seq_length_hist.png".format(args.outDir, name),dpi=600,format="png")
            plt.close()

    # Check that any files were processed
    if len(lengthDict) == 0:
        print("ERROR: No input files were processed...no output generated.")
    else:
        maxLen = 0
        for name in lengthDict:
            lengthDF = lengthDict[name]
            if lengthDF["length"].max() > maxLen:
                maxLen = lengthDF["length"].max()
            plot_dist(lengthDF, "length", name)
        plt.title("Sequence Length Distributions")
        plt.xlim(0, maxLen)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        if args.outPrefix is not None:
            plt.savefig("{}/{}_seq_length_dist.png".format(args.outDir, args.outPrefix),dpi=600,format="png", bbox_inches='tight')
            if args.outPDF:
                plt.savefig("{}/{}_seq_length_dist.pdf".format(args.outDir, args.outPrefix),dpi=600,format="pdf", bbox_inches='tight')
        else:
            plt.savefig("{}/all_seq_length_dist.png".format(args.outDir),dpi=600,format="png", bbox_inches='tight')
            if args.outPDF:
                plt.savefig("{}/all_seq_length_dist.pdf".format(args.outDir),dpi=600,format="pdf", bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
