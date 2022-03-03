#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def restricted_int(val):
    try:
        val = int(val)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not an integer literal" % (val,))
        
    if val < 0:
        raise argparse.ArgumentTypeError("%r not positive integer" % (val,))
    return val

def restricted_float(val):
    try:
        val = float(val)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (val,))

    if val <= 0.0 or val >= 1.0:
        raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0)" % (val,))
    return val

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in a combined rsem expression matrix and output on/off flags based on a minimum value found in at least P% of replicates, or plot values to determine a minimum threshold")

    # Input data
    parser.add_argument("-i", "--input-matrix", dest="inMatrix", required=True, help="Input TSV of combined rsem expression matrix")
    parser.add_argument("-m", "--minimum", dest="min", required=False, type=int, default=0, help="Minimum value required (greater than) in at least P% of replicates to be considered (default: 1)")
    parser.add_argument("-t", "--type", dest="inType", required=False, default='transcript', help="Type of rsem expression file: transcript or gene (default: transcript)")
    parser.add_argument("-u", "--unit", dest="inUnit", required=False, default='TPM', help="Unit of rsem expression expression: TPM or expected_count (default: TPM)")
    parser.add_argument("-e", "--exclude", dest="exclude", required=False, action='append', help="Substring to exclude in sample variables (ex. -e 'noEtoh' to exclude all samples with noEtoh), multiple values can be listed with each '-e'")
    parser.add_argument("-g", "--group", dest="group", required=False, type=restricted_int, action='append', help="Index within sample ID (0_1_2_3...) to group samples by, multiple values can be listed with each '-g' (for example, if sample ID is species_genotype_replicate, -g 1 -g 2 would group all species_genotype replicates, default: 0)")
    parser.add_argument("-p", "--percent", dest="inPercent", type=restricted_float, required=False, default=0.5, help="Percent of replicates in group that are above minimum threshold (default: 0.5, for 50%)")
    parser.add_argument("--plot", dest="plot", required=False, action='store_true', help="Plot distribution of values in rsem expression matrix (x axis will be restricted 500 if max is large)")

    # Output data
    parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Output file prefix for TSV of on/off flags")

    args = parser.parse_args()
    return args

# Function to get unique values of list
def unique(list1):
    x = np.array(list1)
    return np.unique(x)

def main():
    # Get input file type (transcript/gene)
    inType = args.inType
    
    # Get combined expression matrix
    expDF = pd.read_csv(args.inMatrix, sep="\t").set_index(inType+'_id')
    
    # Get minimum value threshold and unit for replicates and
    #     percent of replicates that must meet the threshold
    minVal = args.min
    inUnit = args.inUnit
    percent = args.inPercent
    
    # Get indices for group variables
    group = args.group
    
    # Get output prefix
    outPrefix = args.outPrefix
    
    # Remove excluded columns if provided
    if args.exclude is not None:
        for s in args.exclude:
            expDF = expDF.drop(columns=[c for c in expDF.columns if "_"+s+"_" in c])        
        
    try:
        for flag in unique(["_".join([s.split("_")[i] for i in group]) for s in expDF.columns]):
            values = [f for f in flag.split("_")]
            cols = [col for col in expDF.columns if len(values) == len([value for value in values if value in col])]
            tempDF = expDF[cols].copy()
            numReps = len(list(tempDF))
            # Flag each replicate on/off based on minimum value and percent thresholds specified
            for c in cols:
                tempDF['flag_'+c]= np.where(tempDF[c] > minVal,1,0)
            # Check for at least P% of replicates are flagged to flag group
            tempDF['sum'] = tempDF.filter(regex="flag_").sum(axis=1)
            expDF['flag_'+flag] = np.where(tempDF['sum'] >= numReps*percent,1,0)
            expDF['mean_'+flag] = tempDF[cols].sum(axis=1)/len(cols)
    except IndexError:
        print("ERROR: Index provided is out of range for sample ID variables")
        sys.exit()
        
    # Count transcripts/genes NOT detected in any of the samples
    expDF['sum_flag'] = expDF[[c for c in expDF.columns if "flag_" in c]].sum(axis=1)
    someDF = expDF[expDF['sum_flag']>0].copy()
    print("Frequency of {}s where N number of samples were detected".format(inType))
    print("numSamples\tfreq\n"+expDF['sum_flag'].value_counts(sort=False).to_string())
    print("\nDescriptive values of mean {} for each sample (after removal of {}s not detected in all samples\n{}".format(
            inUnit,inType,someDF[[c for c in someDF.columns if "mean" in c]].describe().to_string()))
    print("\nDescriptive values of mean {} for only detected {}s in each sample:".format(inUnit,inType))
    for column in [c for c in someDF.columns if "mean" in c]:
        print("   ** {} **\n{}".format("_".join(column.split("_")[1:]),someDF[someDF["flag_"+"_".join(column.split("_")[1:])]==1][column].describe().to_string()))

    # Output flags
    expDF.filter(regex='flag_').reset_index().to_csv(outPrefix+"_"+str(minVal)+".tsv", index=False, sep="\t")
    expDF.reset_index().to_csv(outPrefix+"_"+str(minVal)+"_full_table.tsv", index=False, sep="\t")

    # Plot distribution of values
    if args.plot:
        someDF[[c for c in someDF.columns if "mean" in c]].plot.density(figsize=(12,12))
#        plt.xlim(0,round(someDF[[c for c in someDF.columns if "mean" in c]].quantile(0.99).max()))
        plt.xlim(0,min(500,someDF[[c for c in someDF.columns if "mean" in c]].max().max()))
        plt.xlabel("Mean Expression")
        plt.title("Expression values for {0}s detected\nwith > {1} {2} in at least {3:.0%} of replicates".format(inType,minVal,inUnit,percent))
        plt.savefig("{}_{}_{}_{}.png".format(outPrefix,inType,minVal,inUnit),dpi=600,format="png")
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

