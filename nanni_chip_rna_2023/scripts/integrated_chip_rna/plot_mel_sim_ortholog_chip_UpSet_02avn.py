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
            "Plot UpSet plot for the number/proportion of genes that have "
            "H3K4me3/H3K27me2me3 presence for each sex in D. melanoster and D. simulans"
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
    orthoDF1 = pd.read_csv(args.inFile, low_memory=False)

    # Drop orthologs pairs that are not one-2-one
    orthoDF2 = orthoDF1[orthoDF1["flag_one2one_ortholog"]==1]

    # Drop ortholog pairs that the chromosome does not match (e.g. X in mel but A in sim)
    orthoDF = orthoDF2[orthoDF2["xsome"] == orthoDF2["sim_xsome"]].copy()

    # Set up log file if requested
    if args.outLog is not None:
        logfile = open(args.outLog,'w')
        logfile.write("{} total one-to-one orthologous gene pairs".format(
                len(orthoDF))
        )

    # Add column for percent of each gene on X chromosome or Autosomes to be summed in plot
    totalX = len(orthoDF[orthoDF['xsome']=="X"])
    totalA = len(orthoDF[orthoDF['xsome']=="A"])
    orthoDF['Percent'] = np.where(
            orthoDF['xsome']=="X",
            (1.0 / totalX) * 100,
            np.where(
                    orthoDF['xsome']=="A",
                    (1.0 / totalA) * 100,
                    np.nan
            )
    )

    # Convert flags to booleans and give label to use in plot
    flagLabelDict = {
            "flag_any_k4": "Any H3K4me3 Present",
            "flag_has_male_k4": "Male H3K4me3 Present",
            "flag_has_female_k4": "Female H3K4me3 Present",
            "flag_male_limited_k4": "Male-limited H3K4me3",
            "flag_female_limited_k4": "Female-limited H3K4me3",
            "flag_any_k27": "Any H3K27me2me3 Present",
            "flag_has_male_k27": "Male H3K27me2me3 Present",
            "flag_has_female_k27": "Female H3K27me2me3 Present",
            "flag_male_limited_k27": "Male-limited H3K27me2me3",
            "flag_female_limited_k27": "Female-limited H3K27me2me3"
    }
    boolFlags = []
    for species in ["D. melanogaster", "D. simulans"]:
        if species == "D. melanogaster":
            shortName = "mel"
        else:
            shortName = "sim"
        for flag in flagLabelDict.keys():
            orthoDF[species + " " + flagLabelDict[flag]] = orthoDF[flag + "_" + shortName].astype(bool)
            boolFlags.append(species + " " + flagLabelDict[flag])

    # Make list of boolean flags for any H3K4me3 presence flags
    groupListK4 = [c for c in boolFlags if "Any H3K4me3 Present" in c]
    # Make list of boolean flags for just H3K4me3 sex presence flags
    groupListSexK4 = [c for c in boolFlags if "ale H3K4me3 Present" in c]
    # Make list of boolean flags for just male presence H3K4me3 flags
    groupListMaleK4 = [c for c in boolFlags if "Male H3K4me3 Present" in c]
    # Make list of boolean flags for any H3K4me3 presence flags
    groupListK27 = [c for c in boolFlags if "Any H3K27me2me3 Present" in c]
    # Make list of boolean flags for just H3K4me3 sex presence flags
    groupListSexK27 = [c for c in boolFlags if "ale H3K27me2me3 Present" in c]
    # Make list of boolean flags for just male presence H3K4me3 flags
    groupListMaleK27 = [c for c in boolFlags if "Male H3K27me2me3 Present" in c]

    # For each mark
    for chip in ['H3K4me3','H3K27me2me3']:
        if chip == "H3K4me3":
            # Subset for presence flags
            reindexDF = orthoDF.set_index(groupListK4)
            # Subset for sex presence flags
            reindexSexDF = orthoDF.set_index(groupListSexK4)
            # Subset for only those with male-bias in at least one species
            reindexMaleDF = orthoDF.set_index(groupListMaleK4)
        else:
            # Subset for presence flags
            reindexDF = orthoDF.set_index(groupListK27)
            # Subset for sex presence flags
            reindexSexDF = orthoDF.set_index(groupListSexK27)            
            # Subset for only those with male-bias in at least one species
            reindexMaleDF = orthoDF.set_index(groupListMaleK27)
        # For X chromosome or autosomes
        for xsome in ['X','A']:
            # Plot UpSet by degree sort
            # Plot percent fo orthologs on X/A
            upset = UpSet(
                    reindexDF[reindexDF['xsome']==xsome],
                    subset_size='sum',
                    show_counts=True,
                    sum_over='Percent',
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_Percent_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_Percent_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                        dpi=600,
                        format="pdf"
            )
            plt.close()

            # With sex presence flags
            upset = UpSet(
                    reindexSexDF[reindexSexDF['xsome']==xsome],
                    subset_size='sum',
                    show_counts=True,
                    sum_over='Percent',
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_Sex_Percent_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_Sex_Percent_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                        dpi=600,
                        format="pdf"
            )
            plt.close()

            # With male only
            upset = UpSet(
                    reindexMaleDF[reindexMaleDF['xsome']==xsome],
                    subset_size='sum',
                    show_counts=True,
                    sum_over='Percent',
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_hasMale_Percent_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_hasMale_Percent_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="pdf"
            )
            plt.close()

            # Plot count of genes
            upset = UpSet(
                    reindexDF[reindexDF['xsome']==xsome],
                    subset_size='count',
                    show_counts=True,
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_Count_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_Count_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="pdf"
            )
            plt.close()

            # With sex presence flags
            upset = UpSet(
                    reindexSexDF[reindexSexDF['xsome']==xsome],
                    subset_size='count',
                    show_counts=True,
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_Sex_Count_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_Sex_Count_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="pdf"
            )
            plt.close()

            # With male only
            upset = UpSet(
                    reindexMaleDF[reindexMaleDF['xsome']==xsome],
                    subset_size='count',
                    show_counts=True,
                    sort_by='degree',
                    sort_categories_by=None
            )
            upset.plot()
            plt.subplots_adjust(right=1.00001)
            plt.savefig(
                    "{}_{}_UpSet_hasMale_Count_{}.png".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="png"
            )
            plt.savefig(
                    "{}_{}_UpSet_hasMale_Count_{}.pdf".format(
                            args.outPrefix,
                            chip,
                            xsome),
                    dpi=600,
                    format="pdf"
            )
            plt.close()

            # Print percentages and counts to log file if requested
            if args.outLog is not None:
                if chip == "H3K4me3":
                    groupList = groupListK4
                    groupListSex = groupListSexK4
                    groupListMale = groupListMaleK4
                else:
                    groupList = groupListK27
                    groupListSex = groupListSexK27
                    groupListMale = groupListMaleK27

                logfile.write(
                    "Gene Counts ({}, {}):\n{}\n\n{}\n\n{}\n".format(
                        chip,
                        xsome,
                        reindexDF[
                                reindexDF['xsome']==xsome
                            ].reset_index().groupby(groupList)['mel_geneID'].count().to_string(),
                        reindexSexDF[
                                reindexSexDF['xsome']==xsome
                            ].reset_index().groupby(groupListSex)['mel_geneID'].count().to_string(),
                        reindexMaleDF[
                                reindexMaleDF['xsome']==xsome
                            ].reset_index().groupby(groupListMale)['mel_geneID'].count().to_string())
                )
                logfile.write(
                    "\nPercentages ({}, {}):\n{}\n\n{}\n\n{}\n\n".format(
                        chip,
                        xsome,
                        reindexDF[
                                reindexDF['xsome']==xsome
                            ].reset_index().groupby(groupList)['Percent'].sum().to_string(),
                        reindexSexDF[
                                reindexSexDF['xsome']==xsome
                            ].reset_index().groupby(groupListSex)['Percent'].sum().to_string(),
                        reindexMaleDF[
                                reindexMaleDF['xsome']==xsome
                            ].reset_index().groupby(groupListMale)['Percent'].sum().to_string())
                )


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
