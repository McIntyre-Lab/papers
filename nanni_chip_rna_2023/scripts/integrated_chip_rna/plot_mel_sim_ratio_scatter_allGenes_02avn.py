#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Make scatter plot of female/male ratios from exon fragments of all genes "
        "between D. melanogaster and D. simulans."
        )
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="CSV for ratios in all genes with 'mel_ratio' and 'sim_ratio' variables."
    )

    # Output data
    parser.add_argument(
        "-d",
        "--out-directory",
        dest="outDir",
        required=True,
        help="Output directory for plots."
    )

    args = parser.parse_args()
    return args

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

def plot_allGenes(allGeneRatioDf, n_std):
    # Drop conserved male and female biased
    allGeneRatioDf = allGeneRatioDf[
        allGeneRatioDf["flag_conserved_male_bias"] + allGeneRatioDf["flag_conserved_female_bias"]==0
    ]
    # Plot all gene ratios with:
    #   conserved unbiased in gray
    #   female-biased in one species and unbiased in the other in pink
    #   male-biased in one species and unbiased in the other in light blue
    #   male-biased in one species and female-biased in the other in purple
    #   gray x=y line
    fig, ax = plt.subplots(figsize=(14, 7.5))
    sns.scatterplot(
        data=allGeneRatioDf[allGeneRatioDf["flag_F_mel_U_sim"]+allGeneRatioDf["flag_F_sim_U_mel"]==1],
        x="mel_ratio",
        y="sim_ratio",
        alpha=0.3,
        s=20,
        color="red",
        label="Female-biased in one species only",
        ax=ax
    )
    sns.scatterplot(
        data=allGeneRatioDf[allGeneRatioDf["flag_M_mel_U_sim"]+allGeneRatioDf["flag_M_sim_U_mel"]==1],
        x="mel_ratio",
        y="sim_ratio",
        alpha=0.3,
        s=20,
        color="blue",
        label="Male-biased in one species only",
        ax=ax
    )
    sns.scatterplot(
        data=allGeneRatioDf[allGeneRatioDf["flag_conserved_unbiased"]==1],
        x="mel_ratio",
        y="sim_ratio",
        alpha=0.3,
        s=20,
        color="gray",
        label="Conserved unbiased",
        ax=ax
    )
    sns.scatterplot(
        data=allGeneRatioDf[allGeneRatioDf["flag_M_mel_F_sim"]+allGeneRatioDf["flag_M_sim_F_mel"]==1],
        x="mel_ratio",
        y="sim_ratio",
#        alpha=0.3,
        s=20,
        color="black",
#        linewidth=1.5,
#        fc="darkviolet",
        label="Sex reversal (male-biased in one species, female-biased in the other)",
        ax=ax
    )
    plt.xlabel("D. melanogaster ratios")
    plt.ylabel("D. simulans ratios")
    xmin = allGeneRatioDf["mel_ratio"].min()
    ymin = allGeneRatioDf["sim_ratio"].min()
    xmax = allGeneRatioDf["mel_ratio"].max()
    ymax = allGeneRatioDf["sim_ratio"].max()
    ax.plot((0, 0), (min(xmin, ymin),max(xmax, ymax)), color="black")
    ax.plot((min(xmin, ymin),max(xmax, ymax)), (0, 0), color="black")
    ax.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray")
    confidence_ellipse(
        allGeneRatioDf[allGeneRatioDf["flag_M_mel_U_sim"]==1]["mel_ratio"],
        allGeneRatioDf[allGeneRatioDf["flag_M_mel_U_sim"]==1]["sim_ratio"],
        ax,
        n_std=n_std,
        edgecolor='navy',
        lw=3,
        label="95% D. melanogaster male-biased only"
    )
    confidence_ellipse(
        allGeneRatioDf[allGeneRatioDf["flag_M_sim_U_mel"]==1]["mel_ratio"],
        allGeneRatioDf[allGeneRatioDf["flag_M_sim_U_mel"]==1]["sim_ratio"],
        ax,
        n_std=n_std,
        edgecolor='navy',
        lw=3,
        linestyle="--",
        label="95% D. simulans male-biased only"
    )
    confidence_ellipse(
        allGeneRatioDf[allGeneRatioDf["flag_F_mel_U_sim"]==1]["mel_ratio"],
        allGeneRatioDf[allGeneRatioDf["flag_F_mel_U_sim"]==1]["sim_ratio"],
        ax,
        n_std=n_std,
        edgecolor='maroon',
        lw=3,
        label="95% D. melanogaster female-biased only"
    )
    confidence_ellipse(
        allGeneRatioDf[allGeneRatioDf["flag_F_sim_U_mel"]==1]["mel_ratio"],
        allGeneRatioDf[allGeneRatioDf["flag_F_sim_U_mel"]==1]["sim_ratio"],
        ax,
        n_std=n_std,
        edgecolor='maroon',
        lw=3,
        linestyle="--",
        label="95% D. simulans female-biased only"
    )
    ax.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()

    if (len(allGeneRatioDf[allGeneRatioDf["flag_conserved_unbiased"]==1])
        + len(allGeneRatioDf[allGeneRatioDf["flag_F_mel_U_sim"]+allGeneRatioDf["flag_F_sim_U_mel"]==1])
        + len(allGeneRatioDf[allGeneRatioDf["flag_M_mel_U_sim"]+allGeneRatioDf["flag_M_sim_U_mel"]==1])
        + len(allGeneRatioDf[allGeneRatioDf["flag_M_mel_F_sim"]+allGeneRatioDf["flag_M_sim_F_mel"]==1]) != len(allGeneRatioDf)):
        print("WARNING: Genes listed in input file are not as expected.")

def main():
    # Get input file
    allGeneRatioDf = pd.read_csv(args.inFile, low_memory=False)

    # Drop ratios set to 1 or -1
    allGeneRatioDf = allGeneRatioDf[
        (np.abs(allGeneRatioDf["mel_ratio"])!=1)
        & (np.abs(allGeneRatioDf["sim_ratio"])!=1)
    ].copy()

    # Plot all genes
    plot_allGenes(allGeneRatioDf, 2)
    plt.savefig(args.outDir + "/ortho_MF_ratio_all_gene_noConservedMF.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    # Plot all genes excluding genes with consistent mapping bias
    plot_allGenes(allGeneRatioDf[allGeneRatioDf["flag_map_better_2_mel"]+allGeneRatioDf["flag_map_better_2_sim"]==0], 2)
    plt.savefig(args.outDir + "/ortho_MF_ratio_all_gene_noConservedMF_dropConsisMapBias.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
