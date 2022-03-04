#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import sys
import math

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Plot scatter plot given a CSV and two column names."
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help=(
            "Input CSV file columns of values to plot in scatter plot."
        )
    )
    parser.add_argument(
        "-x",
        "--x-column",
        dest="inX",
        required=True,
        help=(
            "Column name in the given CSV for the X-axis values of the scatter plot."
        )
    )
    parser.add_argument(
        "-y",
        "--y-column",
        dest="inY",
        required=True,
        help=(
            "Column name in the given CSV for the Y-axis values of the scatter plot."
        )
    )
    parser.add_argument(
        "-xlab",
        "--x-label",
        dest="inXlab",
        required=False,
        help=(
            "Optional value for X-axis label in the scatter plot."
        )
    )
    parser.add_argument(
        "-ylab",
        "--y-label",
        dest="inYlab",
        required=False,
        help=(
            "Optional value for Y-axis label in the scatter plot."
        )
    )
    parser.add_argument(
        "-xlim",
        "--x-limit",
        dest="inXlim",
        type=float,
        required=False,
        help=(
            "If integer/float values, set limit for X-axis values of the scatter plot."
        )
    )
    parser.add_argument(
        "-ylim",
        "--y-limit",
        dest="inYlim",
        type=float,
        required=False,
        help=(
            "If integer/float values, set limit for Y-axis values of the scatter plot."
        )
    )
    parser.add_argument(
        "-g",
        "--group",
        dest="inGroup",
        required=False,
        help=(
            "Column name in the given CSV for X and Y values to be grouped by "
            "for differenct colors on the scatter plot."
        )
    )
    parser.add_argument(
        "-s",
        "--scale",
        dest="inScale",
        required=False,
        choices=["log10", "log2"],
        help=(
            "Scale the axes by log base 10 (log10) or log base 2 (log2) "
            "(Default with no value: linear scale). NOTE: 1 is added to each "
            "value prior to taking the log to ensure plotting of 0 values."
        )
    )

    # Output data
    parser.add_argument(
        "-r",
        "--pearson-correlation",
        dest="outCorr",
        required=False,
        action= "store_true",
        help="Output Pearson Correlation (r) in bottom right corner of plot (Default: No correlation value)."
    )
    parser.add_argument(
        "-t",
        "--output-type",
        dest="outType",
        choices = ["png", "pdf", "tif", "jpeg", "svg", "ps"],
        required=False,
        default = "png",
        help="Output file type of png, pdf, jepg, svg, ps, or tif (Default: png)."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output scatter plot file name."
    )

    args = parser.parse_args()
    return args

def plot_scatter(df, x, y, xlab=None, ylab=None, xlim=None, ylim=None, group=None, scale=None, correlation=False):
    plt.figure(figsize=(5, 5))
    if scale == "log10":
        df["log10_" + x] = df[x].apply(lambda v: math.log(v+1, 10))
        df["log10_" + y] = df[y].apply(lambda v: math.log(v+1, 10))
        varx = "log10_" + x
        vary = "log10_" + y
    elif scale == "log2":
        df["log2_" + x] = df[x].apply(lambda v: math.log(v+1, 2))
        df["log2_" + y] = df[y].apply(lambda v: math.log(v+1, 2))
        varx = "log2_" + x
        vary = "log2_" + y
    else:
        varx = x
        vary = y
    if group is None:
        sns.scatterplot(
            data=df,
            x=varx,
            y=vary,
            alpha=0.3,
            s=20
        )
    else:
        sns.scatterplot(
            data=df,
            x=x,
            y=y,
            alpha=0.3,
            s=20,
            hue=group
        )
    if xlim is None:
        xlim = df[varx].max()
    plt.xlim(0, xlim)
    if ylim is None:
        ylim = df[vary].max()
    plt.ylim(0, ylim)
    if xlab is None:
        xlab = varx
    plt.xlabel(xlab)
    if ylab is None:
        ylab = vary
    plt.ylabel(ylab)
    if group is not None:
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.plot((0, max(xlim, ylim)), (0, max(xlim, ylim)), color="r")
    if correlation:
        plt.text(
            0.72,
            0.05,
            "r = {0:.4}".format(np.corrcoef(df[varx] ,df[vary])[0][1]),
            transform=plt.gca().transAxes,
            fontsize=14
        )


def main():
    # Set matplotlib parameters to allow for pdf/ps text editing in adobe illustrator
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42

    # Get intput file
    inDF = pd.read_csv(args.inFile, low_memory=False)

    # Check that x and y columns are present
    if args.inX not in inDF.columns:
        print("WARNING: {} column not present in input file.".format(args.inX))
        sys.exit()
    if args.inY not in inDF.columns:
        print("WARNING: {} column not present in input file.".format(args.inY))
        sys.exit()
    if args.inGroup is not None and args.inGroup not in inDF.columns:
        print("WARNING: {} column not present in input file.".format(args.inGroup))
        sys.exit()

    # Plot scatter plot
    plot_scatter(
            inDF,
            args.inX,
            args.inY,
            xlab=args.inXlab,
            ylab=args.inYlab,
            xlim=args.inXlim,
            ylim=args.inYlim,
            group=args.inGroup,
            scale=args.inScale,
            correlation=args.outCorr
    )
    plt.savefig(args.outFile, dpi=600, format=args.outType, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
