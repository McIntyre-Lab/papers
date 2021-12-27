#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Print to stdout counts of genic, intergenic, and border intersections (in that order with spaces between)")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Input TSV file from bedtools intersect with the 5th column having values of genic or intergenic")

    # Output data
    parser.add_argument("--output-counts", dest="outCount", required=True, help="Output file for genic, intergenic, and border interactions (only nonzero counts)")
    parser.add_argument("--output-colors", dest="outColor", required=True, help="Output file for colors that correspond to the counts (for pie chart)")
    parser.add_argument("--output-names", dest="outNames", required=True, help="Output file for names that correspond to the counts (for pie chart)")

    args = parser.parse_args()
    return args

def main():
    
    # Open output files
    outCounts = open(args.outCount,'w')
    outColors = open(args.outColor,'w')
    outNames = open(args.outNames,'w')
    
    # Get intersect file
    bedDF = pd.read_csv(args.inFile,sep="\t",names=['peak_chrom','peak_start','peak_end','peakID','featureType','feature_chrom','feature_start','feature_end','featureID','overlap'])
    
    # Remove duplicate featureType intersections for each peak
    featDF = bedDF[['peakID','featureType']].drop_duplicates()
    
    # Use pivot to set flag values
    featDF['values'] = 1
    collapseDF = featDF.pivot(index='peakID', columns='featureType', values='values').fillna(0)
    
    # Label feature type by flag values
    collapseDF['featureType'] = np.where(collapseDF['genic']+collapseDF['intergenic']==2,"border",
              np.where((collapseDF['genic']==1) & (collapseDF['intergenic']==0),"genic", "intergenic"))

    # Print out counts and colors of genic, intergenic, and border peak intersections (only values that are not zero)
    counts = list()
    colors = list()
    names = list()
    colorDict = {"genic":"C3", "intergenic":"C0", "border":"C4"}
    for feature in ["genic","intergenic","border"]:
        num = len(collapseDF[collapseDF['featureType']==feature])
        if(num>0):
            counts.append(num)
            colors.append(colorDict[feature])
            names.append(feature)
    outCounts.write(",".join(map(str,counts)))
    outColors.write(",".join(colors))
    outNames.write(",".join(names))

    outCounts.close()
    outColors.close()
    outNames.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

