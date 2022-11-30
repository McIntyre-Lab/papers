#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from upsetplot import UpSet

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            "Plot UpSet plot for the number of genes that have male-biased, "
            "female-biased, or unbiased expression in D. melanoster and D. simulans"
        )
    )

    # Input data
    parser.add_argument(
        "-i",
        "--input",
        dest="inFile",
        required=True,
        help="Input CSV file of ortholog combination flags"
    )

    # Output data
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix"
    )
    parser.add_argument(
        "-l",
        "--output-log",
        dest="outLog",
        required=False,
        help="Output log file for counts"
    )

    args = parser.parse_args()
    return args

def main():
    # Get gene inputs
    orthoDF = pd.read_csv(args.inFile, low_memory=False)

    # Drop orthologs pairs that are not one-2-one
    orthoDF1 = orthoDF[orthoDF["flag_one2one_ortholog"]==1]
    
    # Drop ortholog pairs that the chromosome does not match (e.g. X in mel but A in sim)
    orthoDF2 = orthoDF1[orthoDF1["xsome"] == orthoDF1["sim_xsome"]]

    # Drop genes not expressed
    orthoExpress2 = orthoDF2[
            (orthoDF2["flag_expressed_mel"]==1)
            & (orthoDF2["flag_expressed_sim"]==1)
    ].copy()
    
    # Drop sex-limited genes
    orthoExpress = orthoExpress2[
            (orthoExpress2["flag_sex_limited_mel"]==0)
            & (orthoExpress2["flag_sex_limited_sim"]==0)
    ].copy()

    # Set up log file if requested
    if args.outLog is not None:
        logfile = open(args.outLog,'w')
        logfile.write("{} total sex-biased genes (excluding all sex-limited) in both species\n".format(
                len(orthoExpress))
        )
    
    # Add column for percent of each gene on X chromosome or Autosomes to be summed in plot
    totalX = len(orthoExpress[orthoExpress['xsome']=="X"])
    totalA = len(orthoExpress[orthoExpress['xsome']=="A"])
    orthoExpress['Percent'] = np.where(
            orthoExpress['xsome']=="X",
            (1.0 / totalX) * 100,
            np.where(
                    orthoExpress['xsome']=="A",
                    (1.0 / totalA) * 100,
                    np.nan
            )
    )

    # Convert flags to booleans and give label to use in plot
    flagLabelDict = {
            "flag_U": "Unbiased",
            "flag_F": "Female-biased",
            "flag_M": "Male-biased",
            "flag_MF": "Male- and Female-biased",
            "flag_sex_biased": "Sex-biased"
    }
    boolFlags = []
    for species in ["D. melanogaster", "D. simulans"]:
        if species == "D. melanogaster":
            shortName = "mel"
        else:
            shortName = "sim"
        for flag in flagLabelDict.keys():
            orthoExpress[species + " "+flagLabelDict[flag]] = orthoExpress[flag + "_" + shortName].astype(bool)
            boolFlags.append(species + " "+flagLabelDict[flag])

    # Make list of boolean flags for just the male-bias/female-bias/unbias comparison
    # (exclude male and female-bias)
    MFUlist = [c for c in boolFlags if "Sex-biased" not in c and "Male- and Female-biased" not in c]
    MFUlist2 = [c for c in boolFlags if "Sex-biased" not in c]
    # Make list of boolean flags for just sex-biased and unbiased
    sexBlist = [c for c in boolFlags if "Sex-biased" in c or "Unbiased" in c]
    sexBlistSort = [c for c in sexBlist if "Unbiased" in c] + [c for c in sexBlist if "Unbiased" not in c]

    # Plot upset by degree sort
    # Plot number of orthologs on X or A
    upset = UpSet(
            orthoExpress[
                (orthoExpress["xsome"].isin(["X","A"]))
            ].set_index(MFUlist2),
            subset_size='count',
            show_counts=True,
            sort_by='degree',
            sort_categories_by=None
    )
    upset.plot()
    plt.subplots_adjust(right=1.00001)
    plt.savefig("{}_UpSet_Count_XA.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_UpSet_Count_XA.pdf".format(args.outPrefix),dpi=600,format="pdf")
    plt.close()


    # For X chromosome or autosomes
    for xsome in ['X','A']:
        # Plot UpSet by degree sort
        # Plot percent of orthologs on X/A
        # For sex-bias split (excluding male_and_female)
        upset = UpSet(
                orthoExpress[
                    (orthoExpress['xsome']==xsome)
                    & (orthoExpress['flag_MF_mel']==0)
                    & (orthoExpress['flag_MF_sim']==0)
                ].set_index(MFUlist),
                subset_size='sum',
                show_counts=True,
                sum_over='Percent',
                sort_by='degree',
                sort_categories_by=None
        )
        upset.plot()
        plt.subplots_adjust(right=1.00001)
        plt.savefig("{}_UpSet_Percent_{}.png".format(args.outPrefix,xsome),dpi=600,format="png")
        plt.savefig("{}_UpSet_Percent_{}.pdf".format(args.outPrefix,xsome),dpi=600,format="pdf")
        plt.close()

        # For sex-bias combined (including male_and_female) vs. unbiased
        upset = UpSet(
                orthoExpress[
                    (orthoExpress['xsome']==xsome)
                ].set_index(sexBlistSort),
                subset_size='sum',
                show_counts=True,
                sum_over='Percent',
                sort_by='degree',
                sort_categories_by=None
        )
        upset.plot()
        plt.subplots_adjust(right=1.00001)
        plt.savefig("{}_UpSet_sexCombined_Percent_{}.png".format(args.outPrefix,xsome),dpi=600,format="png")
        plt.savefig("{}_UpSet_sexCombined_Percent_{}.pdf".format(args.outPrefix,xsome),dpi=600,format="pdf")
        plt.close()

        # Plot count for genes with sex-bias split (excluding male_and_female)
        upset = UpSet(
                orthoExpress[
                    (orthoExpress['xsome']==xsome)
                    & (orthoExpress['flag_MF_mel']==0)
                    & (orthoExpress['flag_MF_sim']==0)
                ].set_index(MFUlist),
                subset_size='count',
                show_counts=True,
                sort_by='degree',
                sort_categories_by=None
        )
        upset.plot()
        plt.subplots_adjust(right=1.00001)
        plt.savefig("{}_UpSet_Count_{}.png".format(args.outPrefix,xsome),dpi=600,format="png")
        plt.savefig("{}_UpSet_Count_{}.pdf".format(args.outPrefix,xsome),dpi=600,format="pdf")
        plt.close()

        # Plot count for genes with combined sex-bias (including male_and_female)
        upset = UpSet(
                orthoExpress[
                    orthoExpress['xsome']==xsome
                ].set_index(sexBlistSort),
                subset_size='count',
                show_counts=True,
                sort_by='degree',
                sort_categories_by=None
        )
        upset.plot()
        plt.subplots_adjust(right=1.00001)
        plt.savefig("{}_UpSet_sexCombined_Count_{}.png".format(args.outPrefix,xsome),dpi=600,format="png")
        plt.savefig("{}_UpSet_sexCombined_Count_{}.pdf".format(args.outPrefix,xsome),dpi=600,format="pdf")
        plt.close()

        # Print percentages and counts to log file if requested
        if args.outLog is not None:
            logfile.write("Gene Counts ({}):\n{}\n\n{}\n".format(
                xsome,
                orthoExpress[
                        orthoExpress['xsome']==xsome
                    ].groupby(sexBlistSort)['mel_geneID'].count().to_string(),
                orthoExpress[
                        orthoExpress['xsome']==xsome
                    ].groupby(sexBlistSort)['mel_geneID'].count().to_string()
            ))
            logfile.write("\nPercentages ({}):\n{}\n\n{}\n\n".format(
                xsome,
                orthoExpress[
                        orthoExpress['xsome']==xsome
                    ].groupby(sexBlistSort)['Percent'].sum().to_string(),
                orthoExpress[
                        orthoExpress['xsome']==xsome
                    ].groupby(sexBlistSort)['Percent'].sum().to_string()
            ))


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
