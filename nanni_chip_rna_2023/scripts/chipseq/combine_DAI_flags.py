#!/usr/bin/env python3

import pandas as pd
import argparse


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in stacked flags of detection above input in all sample types within a species to make merged flag file")

    # Input data
    parser.add_argument("-i", "--input-stacked", dest="inFlags", required=True, help="Input CSV of flags in stacked format with sampleType and flag_presence columns")

    # Output data
    parser.add_argument("-o", "--outputFile", dest="outFile", required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args
    
def main():
    
    # Get stacked flag file
    stackedDF = pd.read_csv(args.inFlags)
     
    # Group by featureID and pivot flags up by sampleType
    mergedDF = stackedDF.pivot(index='featureID', columns='sampleType', values='flag_presence').reset_index()

    # Output flags
    mergedDF.to_csv(args.outFile, index=False)


   
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

