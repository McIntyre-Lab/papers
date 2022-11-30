#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Make scatter plot of female/male ratios from exon fragments significant "
        "for DE in D. melanogaster vs. D. simulans."
        )
    )

    # Input data
    parser.add_argument(
        "-m",
        "--conserved-male",
        dest="conserveM",
        required=True,
        help="CSV for ratios in conserved male-biased genes."
    )
    parser.add_argument(
        "-f",
        "--conserved-female",
        dest="conserveF",
        required=True,
        help="CSV for ratios in conserved female-biased genes."
    )
    parser.add_argument(
        "--Mmel-Fsim",
        dest="MmelFsim",
        required=True,
        help="CSV for ratios in genes with male-biased expression in D. melanogaster and female-biased in D. simulans."
    )
    parser.add_argument(
        "--Fmel-Msim",
        dest="FmelMsim",
        required=True,
        help="CSV for ratios in genes with female-biased expression in D. melanogaster and male-biased in D. simulans."
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

def plot_scatter_w_chrom(data, x, y, xchrom, ychrom, point_color='red', xlabel="", ylabel="", notPresent = "Not present", both="Both", xonly="X only", yonly="Y only"):
    '''
    Parameters
    ----------
    data : DataFrame
        DataFrame to plot from.
    x : String
        Column name for data to be plotted on X-axis.
    y : String
        Column name for data to be plotted on Y-axis.
    xchrom : String
        Column name for the flag of chromatin presence in the species on the X-axis.
    ychrom : String
        Column name for the flag of chromatin presence in the species on the Y-axis.
    point_color : String, optional
        Color for points on plot without chromatin in either species. The default is 'red'.
    xlabel : String, optional
        X-axis label for plot. The default is "".
    ylabel : String, optional
        Y-axis label for plot. The default is "".
    notPresent : String, optional
        Name in legend for points that do not have chromatin present. The default is "Not present".
    both : String, optional
        Name in legend for points that have chromatin in both species. The default is "Both".
    xOnly : String, optional
        Name in legend for points that have chromatin only in the X species. The default is "X only".
    yOnly : String, optional
        Name in legend for points that have chromatin only in the Y species. The default is "Y only".

    Returns
    -------
    Plot.

    '''
    # Plot scatter plot with:
    #   gray x=y line
    #   points with chromatin mark in both species are purple stars
    #   points with chromatin mark in mel only are yellow squares
    #   points with chromatin mark in sim only are seafoam diamonds
    neitherDf = data[
        data[xchrom] + data[ychrom] == 0
    ]
    bothDf = data[
        data[xchrom] + data[ychrom] == 2
    ]
    xOnlyDf = data[
        (data[xchrom]==1) & (data[ychrom]==0)
    ]
    yOnlyDf = data[
        (data[xchrom]==0) & (data[ychrom]==1)
    ]
    sns.scatterplot(
        data=neitherDf,
        x=x,
        y=y,
        alpha=0.3,
        s=20,
        color=point_color,
        label=notPresent
    )
    sns.scatterplot(
        data=bothDf,
        x=x,
        y=y,
        alpha=0.5,
        s=20,
        color="#6D253E",
        marker="*",
        label=both
    )
    sns.scatterplot(
        data=xOnlyDf,
        x=x,
        y=y,
        alpha=0.5,
        s=20,
        color="#FFC107",
        marker="s",
        label=xonly
    )
    sns.scatterplot(
        data=yOnlyDf,
        x=x,
        y=y,
        alpha=0.5,
        s=20,
        color="#89E3B1",
        marker="d",
        label=yonly
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    xlim = data[x].max()
    ylim = data[y].max()
    plt.plot((0, max(xlim, ylim)), (0, max(xlim, ylim)), color="gray", label="y=x")
    plt.legend()

def plot_regress(data, x, y, point_color='red', ci=0, xlabel="", ylabel="", ttest=False, regColor='blue', xyLine=True):
    '''
    
    Parameters
    ----------
    data : DataFrame
        DataFrame to plot from.
    x : String
        Column name for data to be plotted on X-axis.
    y : String
        Column name for data to be plotted on Y-axis.
    point_color : String, optional
        Color of points outside of confidence interval. The default is 'red'.
    ci : Float, optional
        Confidence interval between 0 and 1 to plot. The default is 0, no CI.
    xlabel : String, optional
        X-axis label for plot. The default is "".
    ylabel : String, optional
        Y-axis label for plot. The default is "".
    ttest : Boolean, optional
        Add t statistic and pvalue to the right of the plot for t-tests of the
        two  nulls: 1) the slope = 0 and 2) the slope = 1. The default is True.
    regColor : String, optional.
        Set color of the regression line. The default is "blue".
    xyLine : Boolean, optional
        Draw gray x=y line. The default is True.

    Returns
    -------
    regress_result : stats.linregress result object.
    Plot.

    '''
    # Plot scatter plot with:
    #   gray x=y line
    #   blue regression line with confidence bands
    #   points outside of confidence bands are red
    # Linear regression
    regress_result = stats.linregress(data[x],data[y])
    if len(data) < 25:
        ci = 0
    # Get critical t value
    criticalT = stats.t.ppf(q=1-(1-ci)/2,df=len(data)-2)
    # Get t statistic given the null that the slope = 1 and two-sided p value
    tval0 = np.abs(regress_result.slope)/regress_result.stderr
    pval0 = stats.t.sf(tval0, len(data)-1)*2
    tval1 = np.abs(regress_result.slope - 1)/regress_result.stderr
    pval1 = stats.t.sf(tval1, len(data)-1)*2
    pvalLT1 = stats.t.sf(tval1, len(data)-1)
    GTCI = data[
        data[y] > (
            regress_result.intercept + (criticalT * regress_result.intercept_stderr) +
            ((regress_result.slope + (criticalT * regress_result.stderr)) * data[x])
        )
    ]
    LTCI = data[
        data[y] < (
            regress_result.intercept - (criticalT * regress_result.intercept_stderr) +
            ((regress_result.slope - (criticalT * regress_result.stderr)) * data[x])
        )
    ]
    EQCI = data[
        (data[y] <= (
            regress_result.intercept + (criticalT * regress_result.intercept_stderr) +
            ((regress_result.slope + (criticalT * regress_result.stderr)) * data[x])
        ))
        & (data[y] >= (
            regress_result.intercept - (criticalT * regress_result.intercept_stderr) +
            ((regress_result.slope - (criticalT * regress_result.stderr)) * data[x])
        ))
    ]
    sns.scatterplot(
        data=GTCI,
        x=x,
        y=y,
        alpha=0.3,
        s=20,
        color=point_color
    )
    sns.scatterplot(
        data=LTCI,
        x=x,
        y=y,
        alpha=0.3,
        s=20,
        color=point_color
    )
    sns.scatterplot(
        data=EQCI,
        x=x,
        y=y,
        alpha=0.3,
        s=20,
        color="gray"
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    xlim = data[x].max()
    ylim = data[y].max()
    if xyLine:
        plt.plot((0, max(xlim, ylim)), (0, max(xlim, ylim)), color="gray", label="y=x")
    if len(data) >= 25:
        plt.plot(
            (0,1),
            (regress_result.intercept,regress_result.intercept+regress_result.slope),
            color=regColor,
            label="y={0:.4f}x+{1:.4f}".format(regress_result.slope,regress_result.intercept)
        )
        if ci != 0:
            plt.fill_between(
                (0,1),
                (regress_result.intercept + (criticalT * regress_result.intercept_stderr), regress_result.intercept + (criticalT * regress_result.intercept_stderr) + regress_result.slope + (criticalT * regress_result.stderr)),
                (regress_result.intercept - (criticalT * regress_result.intercept_stderr), regress_result.intercept - (criticalT * regress_result.intercept_stderr) + regress_result.slope - (criticalT * regress_result.stderr)),
                color="b",
                alpha=0.2,
                label="{}% CI".format(int(ci*100))
            )
        if ttest:
            plt.text(
                1.06,
                0.5,
                "H0: beta1 = 0\n  t={0:.2}\n  pval={1:.2}\n\nH0: beta1 = 1\n  t={2:.2}\n  pval={3:.2}\n\nH0: beta1 >= 1\n  t={2:.2}\n  pval={4:.2}".format(
                    tval0, pval0, tval1, pval1, pvalLT1
            ))
    plt.legend()
    return regress_result

def plot_regress_by_chrom(data, x, y, xname, yname, OSCO="male", point_color='red', ci=0):
    OSCOdict = {"male":
                    {"flag_k4": "flag_has_male_k4_",
                     "name_k4": "male K4",
                     "flag_k27": "flag_has_female_k27_",
                     "name_k27": "female K27"},
                "female":
                    {"flag_k4": "flag_has_female_k4_",
                     "name_k4": "female K4",
                     "flag_k27": "flag_has_male_k27_",
                     "name_k27": "male K27"}
    }
    figMchrom, axMchrom = plt.subplots(nrows=4, ncols=4, figsize=(10,10))
    yNum = 0  # rows
    ## Set x and y limits based on all genes
    xlim = data[x].max()
    ylim = data[y].max()
    for row in axMchrom:
        xNum = 0  # columns
        yNum = yNum + 1
        if yNum == 1:
            subsetPlotY = data[
                (data[OSCOdict[OSCO]["flag_k4"]+yname]==1)
                & (data[OSCOdict[OSCO]["flag_k27"]+yname]==0)
            ]
            yLabel = yname+" "+OSCOdict[OSCO]["name_k4"]
        elif yNum == 2:
            subsetPlotY = data[
                (data[OSCOdict[OSCO]["flag_k4"]+yname]==1)
                & (data[OSCOdict[OSCO]["flag_k27"]+yname]==1)
            ]
            yLabel = yname+" "+OSCOdict[OSCO]["name_k4"]+" "+OSCOdict[OSCO]["name_k27"]
        elif yNum == 3:
            subsetPlotY = data[
                (data[OSCOdict[OSCO]["flag_k4"]+yname]==0)
                & (data[OSCOdict[OSCO]["flag_k27"]+yname]==0)
            ]
            yLabel = yname+" OSCO not present"
        else:
            subsetPlotY = data[
                (data[OSCOdict[OSCO]["flag_k4"]+yname]==0)
                & (data[OSCOdict[OSCO]["flag_k27"]+yname]==1)
            ]
            yLabel = yname+" "+OSCOdict[OSCO]["name_k27"]
        for col in row:
            xNum = xNum + 1
            if xNum == 1:
                subsetPlotXY = subsetPlotY[
                    (subsetPlotY[OSCOdict[OSCO]["flag_k4"]+xname]==1)
                    & (subsetPlotY[OSCOdict[OSCO]["flag_k27"]+xname]==0)
                ]
                xLabel = xname+" "+OSCOdict[OSCO]["name_k4"]
            elif xNum == 2:
                subsetPlotXY = subsetPlotY[
                    (subsetPlotY[OSCOdict[OSCO]["flag_k4"]+xname]==1)
                    & (subsetPlotY[OSCOdict[OSCO]["flag_k27"]+xname]==1)
                ]
                xLabel = xname+" "+OSCOdict[OSCO]["name_k4"]+" "+OSCOdict[OSCO]["name_k27"]
            elif xNum == 3:
                subsetPlotXY = subsetPlotY[
                    (subsetPlotY[OSCOdict[OSCO]["flag_k4"]+xname]==0)
                    & (subsetPlotY[OSCOdict[OSCO]["flag_k27"]+xname]==0)
                ]
                xLabel = xname+" OSCO not present"
            else:
                subsetPlotXY = subsetPlotY[
                    (subsetPlotY[OSCOdict[OSCO]["flag_k4"]+xname]==0)
                    & (subsetPlotY[OSCOdict[OSCO]["flag_k27"]+xname]==1)
                ]
                xLabel = xname+" "+OSCOdict[OSCO]["name_k27"]
            sns.scatterplot(
                data=subsetPlotXY,
                x=x,
                y=y,
                alpha=0.3,
                s=20,
                ax=col,
                color=point_color
            )
            if len(subsetPlotXY) > 25:
                # Linear regression
                regress_result = stats.linregress(subsetPlotXY[x],subsetPlotXY[y])
                col.plot(
                    (0,1),
                    (regress_result.intercept,regress_result.intercept+regress_result.slope),
                    color="b",
                    label="y={0:.4f}x+{1:.4f}".format(regress_result.slope,regress_result.intercept)
                )
                if ci != 0:
                    criticalT = stats.t.ppf(q=1-(1-ci)/2,df=len(subsetPlotXY)-2)
                    col.fill_between(
                        (0,1),
                        (regress_result.intercept + (criticalT * regress_result.intercept_stderr), regress_result.intercept + (criticalT * regress_result.intercept_stderr) + regress_result.slope + (criticalT * regress_result.stderr)),
                        (regress_result.intercept - (criticalT * regress_result.intercept_stderr), regress_result.intercept - (criticalT * regress_result.intercept_stderr) + regress_result.slope - (criticalT * regress_result.stderr)),
                        color="b",
                        alpha=0.2,
                        label="{}% CI".format(int(ci*100))
                    )
                col.legend()
            col.set_xlim(0,xlim)
            col.set_ylim(0,ylim)
            if yNum == 4:
                col.set_xlabel(xLabel)
            else:
                col.set_xlabel("")
            if xNum == 1:
                col.set_ylabel(yLabel)
            else:
                col.set_ylabel("")
            col.plot((0, max(xlim, ylim)), (0, max(xlim, ylim)), color="gray")
    plt.tight_layout()

def main():
    # Get input files
    MconservDf = pd.read_csv(args.conserveM, low_memory=False)
    FconservDf = pd.read_csv(args.conserveF, low_memory=False)
    MmelFsimDf = pd.read_csv(args.MmelFsim, low_memory=False)
    FmelMsimDf = pd.read_csv(args.FmelMsim, low_memory=False)

    # Plot conserved male-biased and conserved female-biased genes with:
    #   gray x=y line
    #   blue regression line with confidence bands
    #   points outside of confidence bands are red
    plot_regress(MconservDf, 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0.95, "D. melanogaster male-biased", "D. simulans male-biased")
    plt.savefig(args.outDir + "/ortho_conserved_M_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(MconservDf, 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0, "D. melanogaster male-biased", "D. simulans male-biased")
    plt.savefig(args.outDir + "/ortho_conserved_M_regress.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(MconservDf[MconservDf['flag_map_better_2_mel']+MconservDf['flag_map_better_2_sim']==0], 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0, "D. melanogaster male-biased", "D. simulans male-biased")
    plt.savefig(args.outDir + "/ortho_conserved_M_regress_dropConsisMapBias.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    plot_regress(FconservDf, 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'red', 0.95, "D. melanogaster female-biased", "D. simulans female-biased", regColor="red")
    plt.savefig(args.outDir + "/ortho_conserved_F_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(FconservDf, 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'red', 0, "D. melanogaster female-biased", "D. simulans female-biased", regColor="red")
    plt.savefig(args.outDir + "/ortho_conserved_F_regress.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(FconservDf[FconservDf['flag_map_better_2_mel']+FconservDf['flag_map_better_2_sim']==0], 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'red', 0, "D. melanogaster female-biased", "D. simulans female-biased")
    plt.savefig(args.outDir + "/ortho_conserved_F_regress_dropConsisMapBias.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    # Plot sex switching genes
    plot_regress(MmelFsimDf, 'one_minus_F_M_ratio_mel', 'one_minus_M_F_ratio_sim', 'blue', 0.95, "D. melanogaster male-biased", "D. simulans female-biased")
    plt.savefig(args.outDir + "/ortho_MmelFsim_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(MmelFsimDf, 'one_minus_F_M_ratio_mel', 'one_minus_M_F_ratio_sim', 'blue', 0, "D. melanogaster male-biased", "D. simulans female-biased")
    plt.savefig(args.outDir + "/ortho_MmelFsim_regress.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    plot_regress(FmelMsimDf, 'one_minus_M_F_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0.95, "D. melanogaster female-biased", "D. simulans male-biased")
    plt.savefig(args.outDir + "/ortho_FmelMsim_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(FmelMsimDf, 'one_minus_M_F_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0, "D. melanogaster female-biased", "D. simulans male-biased")
    plt.savefig(args.outDir + "/ortho_FmelMsim_regress.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()


    # Plot conserved male and female genes on the same plot with ellipses for each
    plot_regress(MconservDf, 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0, "D. melanogaster sex-bias ratio", "D. simulans sex-bias ratio")
    plot_regress(FconservDf, 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'red', 0, "D. melanogaster sex-bias ratio", "D. simulans sex-bias ratio", regColor="red", xyLine=False)
    ax = plt.gca()
    confidence_ellipse(
        MconservDf["one_minus_F_M_ratio_mel"],
        MconservDf["one_minus_F_M_ratio_sim"],
        n_std=2,
        ax=ax,
        edgecolor='navy',
        lw=3
    )
    confidence_ellipse(
        FconservDf["one_minus_M_F_ratio_mel"],
        FconservDf["one_minus_M_F_ratio_sim"],
        n_std=2,
        ax=ax,
        edgecolor='maroon',
        lw=3
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_and_F_regress.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress(MconservDf[MconservDf['flag_map_better_2_mel']+MconservDf['flag_map_better_2_sim']==0], 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'blue', 0, "D. melanogaster sex-bias ratio", "D. simulans sex-bias ratio")
    plot_regress(FconservDf[FconservDf['flag_map_better_2_mel']+FconservDf['flag_map_better_2_sim']==0], 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'red', 0, "D. melanogaster sex-bias ratio", "D. simulans sex-bias ratio", regColor="red", xyLine=False)
    ax = plt.gca()
    confidence_ellipse(
        MconservDf[MconservDf['flag_map_better_2_mel']+MconservDf['flag_map_better_2_sim']==0]["one_minus_F_M_ratio_mel"],
        MconservDf[MconservDf['flag_map_better_2_mel']+MconservDf['flag_map_better_2_sim']==0]["one_minus_F_M_ratio_sim"],
        n_std=2,
        ax=ax,
        edgecolor='navy',
        lw=3
    )
    confidence_ellipse(
        FconservDf[FconservDf['flag_map_better_2_mel']+FconservDf['flag_map_better_2_sim']==0]["one_minus_M_F_ratio_mel"],
        FconservDf[FconservDf['flag_map_better_2_mel']+FconservDf['flag_map_better_2_sim']==0]["one_minus_M_F_ratio_sim"],
        n_std=2,
        ax=ax,
        edgecolor='maroon',
        lw=3
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_and_F_regress_dropConsisMapBias.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    # Plot conserved male-biased and conserved female-biased genes with open chromatin in the given sex for each species
    plot_scatter_w_chrom(
        MconservDf,
        'one_minus_F_M_ratio_mel',
        'one_minus_F_M_ratio_sim',
        'flag_has_male_k4_mel',
        'flag_has_male_k4_sim',
        'blue',
        'D. melanogaster male-biased',
        'D. simulans male-biased',
        'No male K4 present',
        'Male K4 both species',
        'Male K4 Mel only',
        'Male K4 Sim only'
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_w_K4.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_scatter_w_chrom(
        FconservDf,
        'one_minus_M_F_ratio_mel',
        'one_minus_M_F_ratio_sim',
        'flag_has_female_k4_mel',
        'flag_has_female_k4_sim',
        'red',
        'D. melanogaster female-biased',
        'D. simulans female-biased'
        'No female K4 present',
        'Female K4 both species',
        'Female K4 Mel only',
        'Female K4 Sim only'
    )
    plt.savefig(args.outDir + "/ortho_conserved_F_w_K4.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()


    # For sex switching genes with open chromatin in the given sex of the bias
    plot_scatter_w_chrom(
        MmelFsimDf,
        'one_minus_F_M_ratio_mel',
        'one_minus_M_F_ratio_sim',
        'flag_has_male_k4_mel',
        'flag_has_female_k4_sim',
        'blue',
        'D. melanogaster male-biased',
        'D. simulans female-biased'
        'No male (mel) or female (sim) K4 present',
        'Male (mel) and female (sim) K4 present',
        'Male (mel) K4 and no female (sim)',
        'Female (sim) K4 and no male (mel)'
    )
    plt.savefig(args.outDir + "/ortho_MmelFsim_w_K4.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_scatter_w_chrom(
        FmelMsimDf,
        'one_minus_M_F_ratio_mel',
        'one_minus_F_M_ratio_sim',
        'flag_has_female_k4_mel',
        'flag_has_male_k4_sim',
        'blue',
        'D. melanogaster female-biased',
        'D. simulans male-biased'
        'No female (mel) or male (sim) K4 present',
        'Female (mel) and male (sim) K4 present',
        'Female (mel) K4 and no male (sim)',
        'Male (sim) K4 and no female (mel)'
    )
    plt.savefig(args.outDir + "/ortho_FmelMsim_w_K4.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()


    # Plot scatter plot of all genes across the species within orthologous genes
    # Conserved male-biased genes with
    #   1) open chromatin in males, no closed chromatin in females
    #   2) open chromatin in males, closed chromatin in females
    #   3) no open chromatin in males, no closed chromatin in females
    #   4) no open chromatin in males, closed chromatin in females
    # For mel vs sim (16 plots total)
    plot_regress_by_chrom(MconservDf, 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'mel', 'sim', "male", 'blue', 0.95)
    plt.savefig(args.outDir + "/ortho_conserved_M_by_chromatin_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress_by_chrom(MconservDf, 'one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'mel', 'sim', "male", 'blue', 0)
    plt.savefig(args.outDir + "/ortho_conserved_M_by_chromatin.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    plot_regress_by_chrom(FconservDf, 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'mel', 'sim', "female", 'red', 0.95)
    plt.savefig(args.outDir + "/ortho_conserved_F_by_chromatin_wCI.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    plot_regress_by_chrom(FconservDf, 'one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'mel', 'sim', "female", 'red', 0)
    plt.savefig(args.outDir + "/ortho_conserved_F_by_chromatin.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    # Plot scatter plots of genes with k4 in both species and no K27 in either
    MbothK4noK27_regress = plot_regress(
        MconservDf[
            (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
            & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==0)
        ],
        'one_minus_F_M_ratio_mel',
        'one_minus_F_M_ratio_sim',
        'blue',
        0,
        "D. melanogaster male-biased",
        "D. simulans male-biased",
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_maleK4_both_noFemalek27_both.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    FbothK4noK27_regress = plot_regress(
        FconservDf[
            (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
            & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==0)
        ],
        'one_minus_M_F_ratio_mel',
        'one_minus_M_F_ratio_sim',
        'red',
        0,
        "D. melanogaster female-biased",
        "D. simulans female-biased",
        regColor="red"
    )
    plt.savefig(args.outDir + "/ortho_conserved_F_femaleK4_both_noMaleK27_both.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    # Plot scatter plots of genes with k4 in both species and k27 in one species
    MbothK4oneK27_regress = plot_regress(
        MconservDf[
            (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
            & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==1)
        ],
        'one_minus_F_M_ratio_mel',
        'one_minus_F_M_ratio_sim',
        'blue',
        0,
        "D. melanogaster male-biased",
        "D. simulans male-biased",
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_maleK4_both_femaleK27_oneSpecies.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    FbothK4oneK27_regress = plot_regress(
        FconservDf[
            (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
            & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==1)
        ],
        'one_minus_M_F_ratio_mel',
        'one_minus_M_F_ratio_sim',
        'red',
        0,
        "D. melanogaster female-biased",
        "D. simulans female-biased",
        regColor="red"
    )
    plt.savefig(args.outDir + "/ortho_conserved_F_femaleK4_both_maleK27_oneSpecies.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    # Plot scatter plots of genes with k4 in both species and k27 in both species
    MbothK4bothK27_regress = plot_regress(
        MconservDf[
            (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
            & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==2)
        ],
        'one_minus_F_M_ratio_mel',
        'one_minus_F_M_ratio_sim',
        'blue',
        0,
        "D. melanogaster male-biased",
        "D. simulans male-biased",
    )
    plt.savefig(args.outDir + "/ortho_conserved_M_maleK4_both_femaleK27_both.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()
    FbothK4bothK27_regress = plot_regress(
        FconservDf[
            (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
            & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==2)
        ],
        'one_minus_M_F_ratio_mel',
        'one_minus_M_F_ratio_sim',
        'red',
        0,
        "D. melanogaster female-biased",
        "D. simulans female-biased",
        regColor="red"
    )
    plt.savefig(args.outDir + "/ortho_conserved_F_femaleK4_both_maleK27_both.png", dpi=600, format="png", bbox_inches='tight')
    plt.close()

    # Concatenate data for testing regression slopes
    # For male conserved vs female conserved regression, concatenate data into one file
    MconservDf["flag_male"] = 1
    FconservDf["flag_male"] = 0
    MFout = pd.concat([
        MconservDf[['one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'flag_male']].rename(columns={'one_minus_F_M_ratio_mel': 'x', 'one_minus_F_M_ratio_sim': 'y'}),
        FconservDf[['one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'flag_male']].rename(columns={'one_minus_M_F_ratio_mel': 'x', 'one_minus_M_F_ratio_sim': 'y'}),        
    ], ignore_index=False)
    MFout.to_csv(args.outDir + "/ortho_conserved_MF_test_slope.csv", index=False)

    # For male and female with K27 comparisons, make variable cat_chip_k4k27_group
    #   0 = k4 both species, no k27 both species
    #   1 = k4 both species, k27 one species only
    #   2 = k4 both species, k27 both species
    Mconditions = [
        (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
        & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==0),
        (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
        & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==1),
        (MconservDf["flag_has_male_k4_mel"]+MconservDf["flag_has_male_k4_sim"]==2)
        & (MconservDf["flag_has_female_k27_mel"]+MconservDf["flag_has_female_k27_sim"]==2)
    ]
    Mchoices = [0,1,2]
    MconservDf["cat_chip_k4k27_group"] = np.select(Mconditions, Mchoices, "oops")
    Mout = MconservDf[MconservDf["cat_chip_k4k27_group"]!="oops"][['one_minus_F_M_ratio_mel', 'one_minus_F_M_ratio_sim', 'cat_chip_k4k27_group']].rename(columns={'one_minus_F_M_ratio_mel': "x", 'one_minus_F_M_ratio_sim': "y"}).copy()
    Mout.to_csv(args.outDir + "/ortho_conserved_M_K27_test_slope.csv", index=False)
    Fconditions = [
        (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
        & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==0),
        (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
        & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==1),
        (FconservDf["flag_has_female_k4_mel"]+FconservDf["flag_has_female_k4_sim"]==2)
        & (FconservDf["flag_has_male_k27_mel"]+FconservDf["flag_has_male_k27_sim"]==2)
    ]
    Fchoices = [0,1,2]
    FconservDf["cat_chip_k4k27_group"] = np.select(Fconditions, Fchoices, "oops")
    Fout = FconservDf[FconservDf["cat_chip_k4k27_group"]!="oops"][['one_minus_M_F_ratio_mel', 'one_minus_M_F_ratio_sim', 'cat_chip_k4k27_group']].rename(columns={'one_minus_M_F_ratio_mel': "x", 'one_minus_M_F_ratio_sim': "y"}).copy()
    Fout.to_csv(args.outDir + "/ortho_conserved_F_K27_test_slope.csv", index=False)




    # Test differences of regression slopes between:
    #   1) k4 in both and k27 in neither vs. k4 in both and k27 in one species
    #   2) k4 in both and k27 in neither vs. k4 in both and k27 in both species
    # Use Z-test using formula 4 in Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. (1998). Using the correct statistical test for the equality of regression coefficients. Criminology, 36(4), 859–866.
    #   Z = (b1 - b2) / sqrt( (SEb1)^2 + (SEb2)^2 )
    #   (https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1745-9125.1998.tb01268.x#:~:text=A%20frequent%20strategy%20in%20examining%20such%20interactive%20effects,difference%20between%20slopes%20in%20making%20these%20coefficient%20comparisons.)
    # P-value calcuated using two-tailed (Ho: b1 = b2)
    slopeDiff = open(args.outDir + "/ortho_conserved_M_F_K4_K27_slope_diff_tests.txt", "w")

    # 1) k4 in both and k27 in neither vs. k4 in both and k27 in one species
    slopeDiff.write("TEST 1: k4 in both and k27 in neither vs. k4 in both and k27 in one species\n")
    Mz1 = (MbothK4noK27_regress.slope - MbothK4oneK27_regress.slope) / (
            (MbothK4noK27_regress.stderr)**(2) + (MbothK4oneK27_regress.stderr)**(2)
        )**(0.5)
    Mp1 = stats.norm.sf(abs(Mz1))*2
    slopeDiff.write("\tConserved Male - beta1 = {}, beta2 = {}, z = {}, p = {}\n".format(
        MbothK4noK27_regress.slope, MbothK4oneK27_regress.slope, Mz1, Mp1
    ))
    Fz1 = (FbothK4noK27_regress.slope - FbothK4oneK27_regress.slope) / (
            (FbothK4noK27_regress.stderr)**(2) + (FbothK4oneK27_regress.stderr)**(2)
        )**(0.5)
    Fp1 = stats.norm.sf(abs(Fz1))*2
    slopeDiff.write("\tConserved Female - beta1 = {}, beta2 = {}, z = {}, p = {}\n".format(
        FbothK4noK27_regress.slope, FbothK4oneK27_regress.slope, Fz1, Fp1
    ))

    # 2) k4 in both and k27 in neither vs. k4 in both and k27 in both species
    slopeDiff.write("\nTEST 2: k4 in both and k27 in neither vs. k4 in both and k27 in both species\n")
    Mz2 = (MbothK4noK27_regress.slope - MbothK4bothK27_regress.slope) / (
            (MbothK4noK27_regress.stderr)**(2) + (MbothK4bothK27_regress.stderr)**(2)
        )**(0.5)
    Mp2 = stats.norm.sf(abs(Mz2))*2
    slopeDiff.write("\tConserved Male - beta1 = {}, beta2 = {}, z = {}, p = {}\n".format(
        MbothK4noK27_regress.slope, MbothK4bothK27_regress.slope, Mz2, Mp2
    ))
    Fz2 = (FbothK4noK27_regress.slope - FbothK4bothK27_regress.slope) / (
            (FbothK4noK27_regress.stderr)**(2) + (FbothK4bothK27_regress.stderr)**(2)
        )**(0.5)
    Fp2 = stats.norm.sf(abs(Fz2))*2
    slopeDiff.write("\tConserved Female - beta1 = {}, beta2 = {}, z = {}, p = {}\n".format(
        FbothK4noK27_regress.slope, FbothK4bothK27_regress.slope, Fz2, Fp2
    ))
    slopeDiff.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
