#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from upsetplot import UpSet

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Plot UpSet plots of detected and analyzable genes across genotypes."
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input gene flag file."
    )

    # Output data
    parser.add_argument(
        "-d",
        "--drop-none",
        dest="dropNone",
        action="store_true",
        required=False,
        help="Drop genes that are not detected in any genotype."
    )
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Prefix to use for output file and figures."
    )

    args = parser.parse_args()
    return args


def main():
    # Get input file
    flagDF = pd.read_csv(args.inFile, low_memory=False)

    detectFlagCols = []
    for genotype in ["B73", "C123", "Hp301", "Mo17", "NC338"]:
        detectFlagCols.append("flag_detect_"+genotype)
        if len([c for c in flagDF.columns if "flag_detect_ccs" in c or "flag_detect_shrtRd_" in c]) > 0:
            flagDF["flag_detect_"+genotype] = np.where(
                    (flagDF["flag_detect_ccs_"+genotype.lower()+"_gt0"] == 1)
                    | (flagDF["flag_detect_shrtRd_"+genotype.lower()+"_gt0"] == 1),
                    1,
                    0
            )

    # Drop genes not detected in any if requested
    if args.dropNone:
        flagDF = flagDF[flagDF[detectFlagCols].sum(axis=1) > 0]

    # Plot UpSet of detection (detected in ccs or shrtRd in at least one treatment)
    detectCountDF = flagDF[detectFlagCols].copy()
    detectCountDF.columns = detectCountDF.columns.str.split("_").str[2]
    detectConcatDF = pd.concat([flagDF,detectCountDF],axis=1)
    detectConcatDF = detectConcatDF.set_index(list(detectCountDF.columns))
    upset = UpSet(detectConcatDF,subset_size='count',show_counts=True,sort_by='degree')
    upset.plot()
    plt.savefig("{}_detect_atLeast1_trt_upset.png".format(args.outPrefix),dpi=600,format="png")

    # Plot UpSet of analyzable genes (TPM>5 in at least 1 sample of either treatment)
    analyzeCountDF = flagDF[
            [c for c in flagDF.columns if "flag_analyze" in c]
    ].copy()
    analyzeCountDF.columns = analyzeCountDF.columns.str.split("_").str[2].str.title()
    analyzeCountDF = analyzeCountDF.rename(columns={"Nc338": "NC338"})
    analyzeConcatDF = pd.concat([flagDF,analyzeCountDF],axis=1)
    analyzeConcatDF = analyzeConcatDF.set_index(list(analyzeCountDF.columns))
    upset = UpSet(analyzeConcatDF,subset_size='count',show_counts=True,sort_by='degree')
    upset.plot()
    plt.savefig("{}_analyze_atLeast1_trt_upset.png".format(args.outPrefix),dpi=600,format="png")



if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
