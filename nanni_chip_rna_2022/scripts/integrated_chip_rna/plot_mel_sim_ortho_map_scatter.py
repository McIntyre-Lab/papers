#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Make scatter plot of female and male expression in D. melanogaster vs. "
        "D. simulans when both are mapped the same coordinates."
        )
    )

    # Input data
    parser.add_argument(
        "-m",
        "--conserved-male",
        dest="conserveM",
        required=False,
        help="CSV for expression in conserved male-biased genes."
    )
    parser.add_argument(
        "-f",
        "--conserved-female",
        dest="conserveF",
        required=False,
        help="CSV for expression in conserved female-biased genes."
    )
    parser.add_argument(
        "-a",
        "--all",
        dest="all",
        required=False,
        help="CSV for expression in all genes."
    )
    parser.add_argument(
        "-s",
        "--species",
        dest="species",
        required=False,
        choices=["mel", "sim"],
        help=(
            "Species name if the CSV files provided are of a single species "
            "on two sets of coordinates."
        )
    )

    # Output data
    parser.add_argument(
        "-p",
        "--out-prefix",
        dest="outPrefix",
        required=True,
        help="Output directory and prefix for plots."
    )

    args = parser.parse_args()
    return args


def plot_regress(data, x, y, point_color='red', ci=0, xlabel="", ylabel="", ttest=True):
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

    Returns
    -------
    Plot.

    '''
    # Plot scatter plot with:
    #   gray x=y line
    #   blue regression line with confidence bands
    #   points outside of confidence bands are red
    # Linear regression
    regress_result = stats.linregress(data[x],data[y])
    # Get critical t value
    criticalT = stats.t.ppf(q=1-(1-ci)/2,df=len(data)-2)
    # Get t statistic given the null that the slope = 1 and two-sided p value
    tval0 = np.abs(regress_result.slope)/regress_result.stderr
    pval0 = stats.t.sf(tval0, len(data)-1)*2
    tval1 = np.abs(regress_result.slope - 1)/regress_result.stderr
    pval1 = stats.t.sf(tval1, len(data)-1)*2
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
    xmin = data[x].min()
    ymin = data[y].min()
    xmax = data[x].max()
    ymax = data[y].max()
    plt.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray", label="y=x")
    plt.plot(
        (xmin,xmax),
        (regress_result.intercept + (regress_result.slope * xmin),regress_result.intercept + (regress_result.slope * xmax)),
        color="b",
        label="y={0:.4f}x+{1:.4f}".format(regress_result.slope,regress_result.intercept)
    )
    if ci != 0:
        upperB0 = regress_result.intercept + (criticalT * regress_result.intercept_stderr)
        upperB1min = (regress_result.slope + (criticalT * regress_result.stderr)) * xmin
        upperB1max = (regress_result.slope + (criticalT * regress_result.stderr)) * xmax
        lowerB0 = regress_result.intercept - (criticalT * regress_result.intercept_stderr)
        lowerB1min = (regress_result.slope - (criticalT * regress_result.stderr)) * xmin
        lowerB1max = (regress_result.slope - (criticalT * regress_result.stderr)) * xmax
        plt.fill_between(
            (xmin,xmax),
            (upperB0 + upperB1min, upperB0 + upperB1max),
            (lowerB0 + lowerB1min, lowerB0 + lowerB1max),
            color="b",
            alpha=0.2,
            label="{}% CI".format(int(ci*100))
        )
    if ttest:
        plt.text(
            (xmax - xmin) + xmax*0.06,
            (ymax - ymin) * 0.5,
            "H0: beta1 = 0\n  t={0:.2}\n  pval={1:.2}\n\nH0: beta1 = 1\n  t={2:.2}\n  pval={3:.2}".format(
                tval0, pval0, tval1, pval1
        ))
    plt.legend()


def main():
    # Get input files
    if args.conserveM is not None:
        MconservDf = pd.read_csv(args.conserveM, low_memory=False)
    if args.conserveF is not None:
        FconservDf = pd.read_csv(args.conserveF, low_memory=False)
    if args.all is not None:
        allGeneDf = pd.read_csv(args.all, low_memory=False)

    if args.species is not None:
        species = args.species
        if species == "mel":
            fullName = "D. melanogaster"
        else:
            fullName = "D. simulans"
        # Plot conserved male-biased and conserved female-biased genes
        if args.conserveM is not None:
            maleLst = ['avg_m_' + species + "2mel", 'avg_m_' + species + "2sim"]
            plot_regress(MconservDf, maleLst[0], maleLst[1], 'blue', 0.95, fullName+" male-biased (D. melanogaster genome)", fullName+" male-biased (D. simulans genome)")
            plt.savefig(args.outPrefix + "_conserved_M_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(MconservDf, maleLst[0], maleLst[1], 'blue', 0, fullName+" male-biased (D. melanogaster genome)", fullName+" male-biased (D. simulans genome)")
            plt.savefig(args.outPrefix + "_conserved_M_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
    
        if args.conserveF is not None:
            femaleLst = ['avg_f_' + species + "2mel", 'avg_f_' + species + "2sim"]
            plot_regress(FconservDf, femaleLst[0], femaleLst[1], 'red', 0.95, fullName+" female-biased (D. melanogaster genome)", fullName+" female-biased (D. simulans genome)")
            plt.savefig(args.outPrefix + "_conserved_F_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(FconservDf, femaleLst[0], femaleLst[1], 'red', 0, fullName+" female-biased (D. melanogaster genome)", fullName+" female-biased (D. simulans genome)")
            plt.savefig(args.outPrefix + "_conserved_F_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

        if args.all is not None:
            maleLst = ['avg_m_' + species + '2mel', 'avg_m_' + species + '2sim',
                       'flag_map_better_m_' + species + '2mel',
                       'flag_map_better_m_' + species + '2sim']
            femaleLst = ['avg_f_' + species + '2mel', 'avg_f_' + species + '2sim',
                         'flag_map_better_f_' + species + '2mel',
                         'flag_map_better_f_' + species + '2sim']
            ratioLst = ['ratio_' + species + '2mel', 'ratio_' + species + '2sim']
            plot_regress(allGeneDf[maleLst].dropna(), maleLst[0], maleLst[1], 'gray', 0, fullName+" male expression (D. melanogaster genome)", fullName+" male expression (D. simulans genome)")
            plt.savefig(args.outPrefix + "_M_express_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

            if species == "mel":
                colorLst = ['#332288', '#117733']
            else:
                colorLst = ['#117733', '#332288']
            plot_regress(allGeneDf[maleLst].dropna(), maleLst[0], maleLst[1], 'gray', 0, fullName+" male expression (D. melanogaster genome)", fullName+" male expression (D. simulans genome)", ttest=False)
            sns.scatterplot(
                data=allGeneDf[allGeneDf[maleLst[2]]>0][maleLst].dropna(),
                x=maleLst[0],
                y=maleLst[1],
                alpha=0.5,
                s=20,
                color=colorLst[0]
            )
            sns.scatterplot(
                data=allGeneDf[allGeneDf[maleLst[3]]>0][maleLst].dropna(),
                x=maleLst[0],
                y=maleLst[1],
                alpha=0.5,
                s=20,
                color=colorLst[1]
            )
            if species == "mel":
                plt.text(5, -3, len(allGeneDf[allGeneDf[maleLst[2]]>0][maleLst].dropna()), color=colorLst[0])
                plt.text(-5, 10, len(allGeneDf[allGeneDf[maleLst[3]]>0][maleLst].dropna()), color=colorLst[1])
            else:
                plt.text(5, -3, len(allGeneDf[allGeneDf[maleLst[2]]>0][maleLst].dropna()), color=colorLst[0])
                plt.text(-5, 10, len(allGeneDf[allGeneDf[maleLst[3]]>0][maleLst].dropna()), color=colorLst[1])
            plt.xlabel(fullName+" male expression (D. melanogaster genome)")
            plt.ylabel(fullName+" male expression (D. simulans genome)")
            plt.savefig(args.outPrefix + "_M_express_regress_mapBetter2x.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

            plot_regress(allGeneDf[femaleLst].dropna(), femaleLst[0], femaleLst[1], 'gray', 0, fullName+" female expression (D. melanogaster genome)", fullName+" female expression (D. simulans genome)")
            plt.savefig(args.outPrefix + "_F_express_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

            plot_regress(allGeneDf[femaleLst].dropna(), femaleLst[0], femaleLst[1], 'gray', 0, fullName+" female expression (D. melanogaster genome)", fullName+" female expression (D. simulans genome)", ttest=False)
            sns.scatterplot(
                data=allGeneDf[allGeneDf[femaleLst[2]]>0][femaleLst].dropna(),
                x=femaleLst[0],
                y=femaleLst[1],
                alpha=0.5,
                s=20,
                color=colorLst[0]
            )
            sns.scatterplot(
                data=allGeneDf[allGeneDf[femaleLst[3]]>0][femaleLst].dropna(),
                x=femaleLst[0],
                y=femaleLst[1],
                alpha=0.5,
                s=20,
                color=colorLst[1]
            )
            if species == "mel":
                plt.text(5, -3, len(allGeneDf[allGeneDf[femaleLst[2]]>0][femaleLst].dropna()), color=colorLst[0])
                plt.text(-5, 10, len(allGeneDf[allGeneDf[femaleLst[3]]>0][femaleLst].dropna()), color=colorLst[1])
            else:
                plt.text(5, -3, len(allGeneDf[allGeneDf[femaleLst[2]]>0][femaleLst].dropna()), color=colorLst[0])
                plt.text(-5, 10, len(allGeneDf[allGeneDf[femaleLst[3]]>0][femaleLst].dropna()), color=colorLst[1])
            plt.xlabel(fullName+" female expression (D. melanogaster genome)")
            plt.ylabel(fullName+" female expression (D. simulans genome)")
            plt.savefig(args.outPrefix + "_F_express_regress_mapBetter2x.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()   

            # Plot M/F ratio for mapping to each genome
            plot_regress(allGeneDf[ratioLst].dropna(), ratioLst[0], ratioLst[1], 'gray', 0, fullName+" M/F ratio (D. melanogaster genome)", fullName+" M/F ratio (D. simulans genome)", ttest=False)
            plt.savefig(args.outPrefix + "_MF_express_ratio_all.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(allGeneDf[allGeneDf["flag_conserved_male_bias"]==1][ratioLst].dropna(), ratioLst[0], ratioLst[1], 'blue', 0, fullName+" M/F ratio (D. melanogaster genome)", fullName+" M/F ratio (D. simulans genome)", ttest=False)
            plt.savefig(args.outPrefix + "_MF_express_ratio_conserved_M.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(allGeneDf[allGeneDf["flag_conserved_female_bias"]==1][ratioLst].dropna(), ratioLst[0], ratioLst[1], 'red', 0, fullName+" M/F ratio (D. melanogaster genome)", fullName+" M/F ratio (D. simulans genome)", ttest=False)
            plt.savefig(args.outPrefix + "_MF_express_ratio_conserved_F.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            
    else:
        
        # Plot conserved male-biased and conserved female-biased genes with:
        #   gray x=y line
        #   blue regression line with confidence bands
        #   points outside of confidence bands are red
        if args.conserveM is not None:
            plot_regress(MconservDf, 'avg_m_mel', 'avg_m_sim', 'blue', 0.95, "D. melanogaster male-biased", "D. simulans male-biased")
            plt.savefig(args.outPrefix + "_conserved_M_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(MconservDf, 'avg_m_mel', 'avg_m_sim', 'blue', 0, "D. melanogaster male-biased", "D. simulans male-biased")
            plt.savefig(args.outPrefix + "_conserved_M_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

        if args.conserveF is not None:
            plot_regress(FconservDf, 'avg_f_mel', 'avg_f_sim', 'red', 0.95, "D. melanogaster female-biased", "D. simulans female-biased")
            plt.savefig(args.outPrefix + "_conserved_F_regress_wCI.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()
            plot_regress(FconservDf, 'avg_f_mel', 'avg_f_sim', 'red', 0, "D. melanogaster female-biased", "D. simulans female-biased")
            plt.savefig(args.outPrefix + "_conserved_F_regress.png", dpi=600, format="png", bbox_inches='tight')
            plt.close()

    # sns.scatterplot(
    #     data=MconservDf,
    #     x="avg_m_mel",
    #     y="avg_m_sim",
    #     alpha=0.5,
    #     s=20,
    #     color='b'
    # )
    # plt.xlabel("D. melanogaster male-biased")
    # plt.ylabel("D. simulans male-biased")
    # xmin = MconservDf["avg_m_mel"].min()
    # ymin = MconservDf["avg_m_sim"].min()
    # xmax = MconservDf["avg_m_mel"].max()
    # ymax = MconservDf["avg_m_sim"].max()
    # plt.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray", label="y=x")
    

    # sns.scatterplot(
    #     data=FconservDf,
    #     x="avg_f_mel",
    #     y="avg_f_sim",
    #     alpha=0.5,
    #     s=20,
    #     color='r'
    # )
    # plt.xlabel("D. melanogaster female-biased")
    # plt.ylabel("D. simulans female-biased")
    # xmin = FconservDf["avg_f_mel"].min()
    # ymin = FconservDf["avg_f_sim"].min()
    # xmax = FconservDf["avg_f_mel"].max()
    # ymax = FconservDf["avg_f_sim"].max()
    # plt.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray", label="y=x")

    # sns.scatterplot(
    #     data=MmelFsimDf,
    #     x="avg_f_mel",
    #     y="avg_f_sim",
    #     alpha=0.5,
    #     s=20,
    # )
    # plt.xlabel("D. melanogaster male-biased")
    # plt.ylabel("D. simulans female-biased")
    # xmin = MmelFsimDf["avg_m_mel"].min()
    # ymin = MmelFsimDf["avg_m_sim"].min()
    # xmax = MmelFsimDf["avg_m_mel"].max()
    # ymax = MmelFsimDf["avg_m_sim"].max()
    # plt.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray", label="y=x")

    # sns.scatterplot(
    #     data=FmelMsimDf,
    #     x="avg_f_mel",
    #     y="avg_m_sim",
    #     alpha=0.5,
    #     s=20,
    # )
    # plt.xlabel("D. melanogaster female-biased")
    # plt.ylabel("D. simulans male-biased")
    # xmin = FmelMsimDf["avg_f_mel"].min()
    # ymin = FmelMsimDf["avg_m_sim"].min()
    # xmax = FmelMsimDf["avg_f_mel"].max()
    # ymax = FmelMsimDf["avg_m_sim"].max()
    # plt.plot((min(xmin, ymin), max(xmax, ymax)), (min(xmin, ymin), max(xmax, ymax)), color="gray", label="y=x")




if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
