#!/usr/bin/env python3

import pandas as pd
import numpy as np
import itertools
import argparse


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in FRiP values across feature types for all sample IDs within a species and creates table of feature FRiP ratios")

    # Input data
    parser.add_argument("-i", "--input-frip", dest="inFRiP", required=True, help="Input CSV of FRiP values with all sampleID and featureTypes")

    # Output data
    parser.add_argument("-o", "--outputFile", dest="outFile", required=True, help="Output CSV file name for CSV of all FRiP values per sampleType and proportions of each")
    parser.add_argument("-s", "--outputStat", dest="outStat", required=True, help="Output CSV of min, median, and max values per antibody for each feature ratio")

    args = parser.parse_args()
    return args
    
def main():
    
    # Get frip file
    fripDF = pd.read_csv(args.inFRiP)
     
    # Group by sampleID and pivot frip up by featureType
    sampleDF = fripDF.pivot(index='sampleID', columns='feature_type', values='frip')
    comboDF = pd.DataFrame({'{}_{}'.format(first,second): sampleDF[first]/sampleDF[second] for first,second in itertools.combinations(sampleDF,2)})

    # Output ratios
    comboDF.to_csv(args.outFile)
    
    # Add antibody column
    comboDF['antibody'] = np.where(comboDF.index.str.contains("K4"),"K4",np.where(comboDF.index.str.contains("K27"),"K27","I"))
    
    # Output min, median, and max for each antibody to check for major outliers
    abDF = comboDF.groupby('antibody').agg(['min','median','max']).unstack(1).reset_index().rename(columns={'level_0':'ratio','level_1':'stat',0:'value'}).sort_values(by=['antibody','ratio'])
    abDF.to_csv(args.outStat, index=False)

   
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

