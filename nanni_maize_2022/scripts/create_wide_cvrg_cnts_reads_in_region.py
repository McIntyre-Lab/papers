#!/usr/bin/env python

import argparse
import os
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Sum coverage counts reads in region.")

    # Input data
    parser.add_argument(
        "-d",
        "--directory",
        dest="inD",
        required=True,
        help="Input directory for containing cvrg_cnts_gene_* files to be summed."
    )
    parser.add_argument(
        "-n",
        "--name",
        dest="inN",
        required=True,
        help="Variable name to merge files on."
    )
    parser.add_argument(
        "-g",
        "--genotypes",
        dest="inG",
        required=True,
        help="Comma separated list of genotypes to merge (e.g. B73,Mo17,NC338)."
    )

    # Output data
    parser.add_argument(
        "-p",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix."
    )

    args = parser.parse_args()
    return args

def list_files(directory, suffix, allow_ext=False, prefix=None):
    if prefix is None:
        if allow_ext:
            return (f for f in os.listdir(directory) if f.endswith(suffix) or suffix+"." in f)
        else:
            return (f for f in os.listdir(directory) if f.endswith(suffix))
    else:
        if allow_ext:
            return (f for f in os.listdir(directory) if (f.endswith(suffix) or suffix+"." in f) and (f.startswith(prefix)))
        else:
            return (f for f in os.listdir(directory) if (f.endswith(suffix)) and (f.startswith(prefix)))

def main():
    inD = args.inD
    mergeVar = args.inN
    genotypes = args.inG.split(",")
    # Merge all files in input directory by gene for each genotype given
    for genotype in genotypes:
        first = True
        for file in list_files(inD, ".csv", prefix="cvrg_cnts_gene_"+genotype):
            name = "_".join(file[:-4].split("_")[3:7])
            if "AMB" in name:
                name = name[:-2] + "mb"
            if first:
                fullDF = pd.read_csv(inD+"/"+file, low_memory=False)
                fullDF.columns = fullDF.columns + "_" + name
                fullDF = fullDF.rename(columns={mergeVar+"_"+name: mergeVar})
                first = False
            else:
                tempDF = pd.read_csv(inD+"/"+file, low_memory=False)
                tempDF.columns = tempDF.columns + "_" + name
                tempDF = tempDF.rename(columns={mergeVar+"_"+name: mergeVar})
                tempFullDF = pd.merge(
                    fullDF,
                    tempDF,
                    how="outer",
                    on=mergeVar,
                    validate="1:1"
                )
                fullDF = tempFullDF.copy()
        # Sum mapped reads across replicates of the same treatment within each genotype
        fullDF["sum_reads_in_region_"+genotype+"_Amb"] = fullDF[
                [c for c in fullDF.columns if "reads_in_region" in c and "_Amb" in c]
            ].sum(axis=1)
        fullDF["sum_reads_in_region_"+genotype+"_Ele"] = fullDF[
                [c for c in fullDF.columns if "reads_in_region" in c and "_Ele" in c]
            ].sum(axis=1)
        # Output mapped reads for each replicate and the sum for the genotype
        readCols = [c for c in fullDF.columns if "reads_in_region" in c]
        fullDF[[mergeVar]+readCols].to_csv(
                "{}_reads_in_region_{}_sbys.csv".format(args.outPrefix, genotype),
                index=False
        )

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
