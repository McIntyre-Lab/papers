#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Plot gene group distributions of values.")

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input Tajima's D or H12 for each gene (must have gene_id column)."
    )
    parser.add_argument(
        "-v",
        "--variable",
        dest="inVar",
        required=True,
        help="Input variable name for value to be plot (e.g., tajimaD or H12)."
    )
    parser.add_argument(
        "-n",
        "--name",
        dest="inName",
        required=True,
        help="Input variable name to use in plot for value to be plotting (e.g., Tajima's D or H12')."
    )
    parser.add_argument(
        "-g",
        "--gene-flags",
        dest="inFlag",
        required=True,
        help="Input CSV of orthologous gene results."
    )
    

    # Output data
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output directory and prefix for plots."
    )

    args = parser.parse_args()
    return args

def main():
    # Get input file and species name
    dataDf = pd.read_csv(args.inFile, low_memory=False)
    geneDf = pd.read_csv(args.inFlag, low_memory=False)
    inVar = args.inVar
    inName = args.inName

    # Get groups of one-to-one genes of interest:
    #   conserved male-biased
    #   conserved femelae-biased
    #   conserved unbiased
    #   male-biased in one species and unbiased in the other
    #   female-biased in one species and unbiased in the other
    #   male-biased in one species and female-biased in the other
    one2oneDf = geneDf[geneDf["flag_one2one_XA_match"]==1].copy()
    geneGroupConditions = [
        one2oneDf["flag_conserved_male_bias"] ==1,
        one2oneDf["flag_conserved_female_bias"]==1,
        one2oneDf["flag_conserved_unbiased"]==1,
        one2oneDf["flag_M_mel_U_sim"]
        + one2oneDf["flag_M_sim_U_mel"]
        + one2oneDf["flag_F_mel_U_sim"]
        + one2oneDf["flag_F_sim_U_mel"]
        + one2oneDf["flag_M_mel_F_sim"]
        + one2oneDf["flag_M_sim_F_mel"] > 0,
    ]
    geneGroupChoices = [
        "Conserved male-biased",
        "Conserved female-biased",
        "Conserved unbiased",
        "Sex switching"
    ]
    one2oneDf["gene_group"] = np.select(geneGroupConditions, geneGroupChoices, "Other")

    # Merge in tajima's D values
    geneTD = pd.merge(
        one2oneDf,
        dataDf,
        how='outer',
        left_on='sim_geneID',
        right_on='gene_id',
        indicator='merge_check',
        validate='1:m'
    )
    geneTD = geneTD[geneTD["merge_check"]=="both"].copy()
    geneTD = geneTD[geneTD["gene_group"]!="Other"].copy()
    geneTD["gene_group"] = pd.Categorical(geneTD["gene_group"], geneGroupChoices)
    geneTD = geneTD.sort_values('gene_group')
    geneTD[inVar] = geneTD[inVar].astype(float)

    # Make density plot of tajima's D for each group for X and A combined
#    plt.figure(figsize=(8,8))
#    ax = sns.distplot(geneTD, x=inVar, hue=geneTD['gene_group'], kind='kde')
#    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
#    plt.xlabel(inName)
#    plt.savefig("{}_gene_{}_density.png".format(args.outPrefix, inVar), dpi=600, format="png", bbox_inches='tight')
#    plt.close()    

    # Plot boxplot of tajima's D for each group for X and A separately and combined
    boxplotData = []
    labelLst = []
    plt.figure(figsize=(10,5))
    for chrom in ["X", "A"]:
        for group in geneTD["gene_group"].unique():
            boxplotData.append(geneTD[(geneTD["gene_group"]==group)&(geneTD["xsome"]==chrom)][inVar])
            labelLst.append(group+" "+chrom)
    plt.boxplot(boxplotData, labels=labelLst, patch_artist=True, vert=False)
    plt.tight_layout()
    plt.savefig("{}_gene_{}_boxplot_by_chrom.png".format(args.outPrefix, inVar), dpi=600, format="png", bbox_inches='tight')
    plt.close()

    boxplotData = []
    labelLst = []
    plt.figure(figsize=(10,5))
    for group in geneTD["gene_group"].unique():
        boxplotData.append(geneTD[(geneTD["gene_group"]==group)][inVar])
        labelLst.append(group)
    plt.boxplot(boxplotData, labels=labelLst, patch_artist=True, vert=False)
    plt.tight_layout()
    plt.savefig("{}_gene_{}_boxplot.png".format(args.outPrefix, inVar), dpi=600, format="png", bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
