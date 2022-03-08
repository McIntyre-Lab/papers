#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Evaluate and plot transcript lengths from FASTA."
    )

    # Input data
    parser.add_argument(
        "-fa",
        "--fasta",
        dest="inFA",
        required=True,
        help="Input FA file of interest."
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix."
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


def main():
    # Get intput FASTA
    fasta = args.inFA

    # Read in FA sequences and store lengths
    lengthDF = pd.DataFrame(columns=["name", "length"])
    for record in SeqIO.parse(fasta, "fasta"):
        lengthDF = lengthDF.append(
            pd.DataFrame(
                    [[record.id, len(record.seq)]], columns=["name", "length"]
            ),
            ignore_index=True
        )
    lengthDF.to_csv(args.outPrefix + "_seq_length.csv", index=False)

    # Plot sequence length distribution
    plot_hist(lengthDF, "length", "Sequence Length Distribution")
    plt.savefig(args.outPrefix + "_seq_length_hist.png",dpi=600,format="png")
    plt.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
