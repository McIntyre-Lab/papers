#!/usr/bin/env python

import argparse
import pandas as pd
from matplotlib import pyplot as plt

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Plot bar plots of expression and chromatin relationship."
        )
    )

    # Input data
    parser.add_argument(
        "-m",
        "--mel",
        dest="inMG",
        required=True,
        help=(
            "Input CSV of D. melanogaster gene-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-s",
        "--sim",
        dest="inSG",
        required=True,
        help=(
            "Input CSV of D. simulans gene-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-o",
        "--ortho",
        dest="inO",
        required=True,
        help=(
            "Input CSV of D. melangoster and D. simulans orthologous "
            "gene-level file (from supplement)."
        )
    )

    # Output data
    parser.add_argument(
        "-d",
        "--output-directory",
        dest="outDir",
        required=True,
        help="Output direcotry"
    )

    args = parser.parse_args()
    return args

def main():
    # Get gene inputs
    melDF = pd.read_csv(args.inMG, low_memory=False)
    simDF = pd.read_csv(args.inSG, low_memory=False)
    orthoDF = pd.read_csv(args.inO, low_memory=False)

    # Drop orthologs pairs that are not one-2-one
    orthoDF = orthoDF[orthoDF["flag_one2one_XA_match"]==1]
    melOrthoGenes = orthoDF["mel_geneID"].unique()
    simOrthoGenes = orthoDF["sim_geneID"].unique()

    # Drop genes not expressed
    melExpress2 = melDF[melDF["flag_expressed"] == 1].copy()
    simExpress2 = simDF[simDF["flag_expressed"] == 1].copy()
    
    # Drop sex-limited genes
    melExpress = melExpress2[melExpress2["flag_sex_limited"] == 0]
    simExpress = simExpress2[simExpress2["flag_sex_limited"] == 0]

    # Get dataframe subsets of X and A (no chrom 4) only
    melExpressXA = melExpress[melExpress["xsome"].isin(["X", "A"])]
    simExpressXA = simExpress[simExpress["xsome"].isin(["X", "A"])]

    # Loop over all genes vs. subset orthologs
    for group in ["all", "orthologous"]:
        # Loop over species
        for species in ["mel", "sim"]:
            if species == "mel":
                if group == "orthologous":
                    df = melExpressXA[melExpressXA["FBgn"].isin(melOrthoGenes)]
                else:
                    df = melExpressXA
            else:
                if group == "orthologous":
                    df = simExpressXA[simExpressXA["fbgn"].isin(simOrthoGenes)]
                else:
                    df = simExpressXA
            # Loop over X and autosome
            for xsome in ["X", "A", "XA"]:
                if xsome == "XA":
                    chromExpress = df
                else:
                    chromExpress = df[df["xsome"]==xsome]
                # Loop over sex-biased expression classification
                for bias in ["Male-biased", "Female-biased"]:
                    # Get associated expression and chromatin flags
                    if bias == "Male-biased":
                        expFlag = "flag_M"
                        openFlag = "flag_has_male_k4"
                        closeFlag = "flag_has_female_k27"
                        color = "b"
                    else:
                        expFlag = "flag_F"
                        openFlag = "flag_has_female_k4"
                        closeFlag = "flag_has_male_k27"
                        color = "r"
                    # Plot percent of sex-biased genes with open chromatin present
                    barwidth = 0.8
                    # plt.figure(figsize=(1.75,2))
                    plt.figure(figsize=(1.5,2))
                    biasOpen = (
                            chromExpress[chromExpress[expFlag]==1][openFlag].sum() / chromExpress[expFlag].sum()
                        ) * 100
                    plt.bar(1,
                            biasOpen,
                            color=color,
                            edgecolor='w',
                            width = barwidth
                    )
                    plt.text(1, biasOpen + 1, round(biasOpen, 2), ha="center")
                    # Plot percent of gene without given sex-bias with open chromatin present
                    noBiasOpen = (
                            chromExpress[chromExpress[expFlag]!=1][openFlag].sum() / len(chromExpress[chromExpress[expFlag]!=1])
                        ) * 100
                    plt.bar(2,
                            noBiasOpen,
                            color=color,
                            edgecolor='w',
                            hatch="/",
                            # color='k',
                            width = barwidth
                    )
                    plt.text(2, noBiasOpen + 1, round(noBiasOpen, 2), ha="center")
                    plt.xticks([])
                    plt.ylim(0,100)
                    # plt.ylabel("% Genes with (solid) or without (hatched) \n{} Expression".format(bias))
                    # plt.ylabel("% Genes")
                    plt.tight_layout()
                    plt.savefig("{}/{}_{}_{}_{}_open.png".format(args.outDir, group, species, xsome, bias), dpi=600, format="png")
                    plt.close()
    
                    
                    # Plot percent of sex-biased genes with closed chromatin present
                    # plt.figure(figsize=(1.75,2))
                    plt.figure(figsize=(1.5,2))
                    biasClosed = (
                            chromExpress[chromExpress[expFlag]==1][closeFlag].sum() / chromExpress[expFlag].sum()
                        ) * 100
                    plt.bar(1,
                            biasClosed,
                            color=color,
                            edgecolor='w',
                            width = barwidth
                    )
                    plt.text(1, biasClosed + 1, round(biasClosed, 2), ha="center")
                    # Plot percent of gene without given sex-bias with closed chromatin present
                    noBiasClosed = (
                            chromExpress[chromExpress[expFlag]!=1][closeFlag].sum() / len(chromExpress[chromExpress[expFlag]!=1])
                        ) * 100
                    plt.bar(2,
                            noBiasClosed,
                            color=color,
                            edgecolor='w',
                            hatch="/",
                            # color='k',
                            width=barwidth
                    )
                    plt.text(2, noBiasClosed + 1, round(noBiasClosed, 2), ha="center")
                    plt.xticks([])
                    plt.ylim(0,100)
                    # plt.ylabel("% Genes with (solid) or without (hatched) \n{} Expression".format(bias))
                    # plt.ylabel("% Genes")
                    plt.tight_layout()
                    plt.savefig("{}/{}_{}_{}_{}_closed.png".format(args.outDir, group, species, xsome, bias), dpi=600, format="png")
                    plt.close()


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
