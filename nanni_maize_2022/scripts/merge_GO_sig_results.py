#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Merge significant GO term files by GO term."
    )

    # Input data
    parser.add_argument(
        "-f",
        "--files",
        dest="inFiles",
        required=True,
        help=(
            "Comma separated list of full file paths to file with significant "
            "GO terms (all files must contain 'Term' column)."
        )
    )
    parser.add_argument(
        "-n",
        "--names",
        dest="inNames",
        required=True,
        help=(
            "Comma separated list of names corresponding to the files given in "
            "-f. These names will be used as prefixes to columns and must be "
            "in the same order as the files."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output",
        dest="outFile",
        required=True,
        help="Output CSV file for combined GO results."
    )
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=False,
        help="Output prefix for compisons of GO results."
    )
    parser.add_argument(
        "-c",
        "--comparison",
        dest="outComp",
        required=False,
        action="append",
        help=(
            "Comma separated pair of names in -n to compare GO results of. "
            "Each additional -c argument will add a new comparison - "
            "Must have -p option for outputs or no comparisons will be made."
        )
    )

    args = parser.parse_args()
    return args

def main():
    # Splt lists of files and names
    files = args.inFiles.split(",")
    names = args.inNames.split(",")
    if len(files) != len(names):
        print("!!!ERROR: Different number of files and names.")
        exit()
    first = True
    for ind in range(0,len(files)):
        if first:
            # Get file and add prefix to all columns
            combGO = pd.read_csv(files[ind])
            # Update columns with flag_detect_DE to be flag_DE and all5 to be all
            for col in combGO.columns:
                if "flag_detect" in col:
                    combGO = combGO.rename(columns={
                            col: "flag_" + col.split("flag_detect_")[1]
                    })
                    col = "flag_" + col.split("flag_detect_")[1]
                if "all5" in col:
                    combGO = combGO.rename(columns={
                            col: col.split("all5")[0] + "all"
                    })
            # Verify terms were read in properly and all flag columns as 0 or 1
            flagCols = [c for c in combGO.columns if c != "Term"]
            if not combGO[flagCols].isin([0,1]).all().all():
                print(
                    "!!!ERROR: Input file not read in properly and contains "
                    "flag column that is not 0/1:\n{}\n".format(
                            files[ind]
                    )
                )
            combGO.columns = names[ind] + "_" + combGO.columns
            combGO = combGO.rename(columns={names[ind] + "_Term": "Term"})
            first = False
        else:
            # Get file and add prefix to all columns
            tempGO = pd.read_csv(files[ind])
            # Update columns with flag_detect_DE to be flag_DE
            for col in combGO.columns:
                if "flag_detect" in col:
                    combGO = combGO.rename(columns={
                            col: "flag_" + col.split("flag_detect_")[1]
                    })
                    col = "flag_" + col.split("flag_detect_")[1]
                if "all5" in col:
                    combGO = combGO.rename(columns={
                            col: col.split("all5")[0] + "all"
                    })
            # Verify terms were read in properly and all flag columns as 0 or 1
            flagCols = [c for c in tempGO.columns if c != "Term"]
            if not tempGO[flagCols].isin([0,1]).all().all():
                print(
                    "!!!ERROR: Input file not read in properly and contains "
                    "flag column that is not 0/1:\n{}\n".format(
                            files[ind]
                    )
                )
            tempGO.columns = names[ind] + "_" + tempGO.columns
            tempGO = tempGO.rename(columns={names[ind] + "_Term": "Term"})

            # Merge with combined file
            tempMergeGO = pd.merge(
                    combGO,
                    tempGO,
                    how="outer",
                    on="Term",
                    validate="1:1",
            )
            combGO = tempMergeGO.fillna(0).copy()

    # Evaluate similarities and differences between each enrichment if requested
    if args.outPrefix is not None:
        # Count all similarites for each flag
        similarCountDF = pd.DataFrame()
        for col in flagCols:
            groupCols = [c for c in combGO.columns if c.endswith(col)]
            if len(groupCols) < len(files):
                print("!!!WARNING: {} in {} of the {} enrichment tests".format(
                        col,
                        len(groupCols),
                        len(files)
                ))
            groupTermDF = combGO[["Term"]].copy()
            groupTermDF[col] = combGO[groupCols].sum(axis=1).copy()
            if groupTermDF[col].sum() > 0:
                if len(similarCountDF) == 0:
                    similarCountDF = groupTermDF.copy()
                else:
                    tempCountMerge = pd.merge(
                            similarCountDF,
                            groupTermDF[groupTermDF[col]>0],
                            how="outer",
                            on="Term",
                            validate="1:1"
                    )
                    similarCountDF = tempCountMerge.copy().fillna(0)
        similarCountDF.to_csv("{}_identical_all_{}_enrichments.csv".format(
                    args.outPrefix,
                    len(files)),
                index=False
        )
        # Do comparisons requested
        if args.outComp is not None:
            for comp in args.outComp:
                if len(comp.split(",")) != 2:
                    print(
                        "!!!WARNING: More than 2 names provided for comparison "
                        "{} - only using first 2".format(comp)
                    )
                name1 = comp.split(",")[0]
                name2 = comp.split(",")[1]
                if name1 not in names or name2 not in names:
                    print(
                        "!!!ERROR: Comparison names not recognized in {} "
                        "- skipping.".format(
                                comp
                        )
                    )
                    continue
                compareDF = pd.DataFrame()
                for col in flagCols:
                    groupCols = [name1 + "_" + col, name2 + "_" + col]
                    
                    if name1 + "_" + col not in combGO.columns:
                        print("!!!WARNING: {} only not in {}...skipping.".format(
                                col,
                                name1
                        ))
                        continue
                    if name2 + "_" + col not in combGO.columns:
                        print("!!!WARNING: {} only not in {}...skipping.".format(
                                col,
                                name2
                        ))
                        continue
                    groupTermDF = combGO[["Term"]].copy()
                    groupConditions = [
                        combGO[groupCols].sum(axis=1) == 2,
                        (combGO[groupCols].sum(axis=1) == 1)
                        & (combGO[name1 + "_" + col] == 1),
                        (combGO[groupCols].sum(axis=1) == 1)
                        & (combGO[name2 + "_" + col] == 1)
                    ]
                    groupValues = ["2", name1, name2]
                    groupTermDF[col] = np.select(groupConditions, groupValues, "0")
                    if len(groupTermDF[groupTermDF[col]!="0"]) > 0:
                        if len(compareDF) == 0:
                            compareDF = groupTermDF.copy()
                        else:
                            tempCompareMerge = pd.merge(
                                    compareDF,
                                    groupTermDF[groupTermDF[col]!="0"],
                                    how="outer",
                                    on="Term",
                                    validate="1:1"
                            )
                            compareDF = tempCompareMerge.copy().fillna(0)
                compareDF.to_csv("{}_identical_{}_{}_enrichments.csv".format(
                            args.outPrefix,
                            name1,
                            name2),
                        index=False
                )

    # Drop any columns that are all 0 (no significant GO term)
    # Output final combined file
    combGO.loc[:, (combGO != 0).any(axis=0)].to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
