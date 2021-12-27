#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

## Take in the stacked coverage count file of each sample type and 
##     calculate FRiP from reads_in_region and mapped_reads columns
## FRiP is fraction of reads in peaks developed by Ji et al., 2008
##   (Num_reads_in_peaks)/(Num_mapped_reads) * 100

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in a single stacked counts matrix (cat multiple together) to calculate Fraction of Reads in Peaks (FRiP) or (Num_reads_in_peaks)/(Num_mapped_reads) * 100")

    # Input data
    parser.add_argument("-i", "--input-stacked", dest="inCounts", required=True, help="Input CSV of counts in stacked format")

    # Output data
    parser.add_argument("-o", "--output", dest="outTSV", required=True, help="Output file name for TSV of FRiP values")

    args = parser.parse_args()
    return args
def main():  
    
    # Import counts data
    countsStackedDF = pd.read_csv(args.inCounts)
    countsStackedDF = countsStackedDF.rename(columns={'fusion_id':'peakID'})
    print("Samples in Counts File:\n\n"+','.join(countsStackedDF['sampleID'].unique())+"\n\n")
    
    # Check that each sampleID has the same number of peaks
    count = len(countsStackedDF.groupby(['sampleID'])['apn'].count().unique())
    if count != 1 :
        print("ERROR : Number of Peaks is not the same for all samples")
        sys.exit()
    
    # Get mapped_reads for each sample and sum all reads_in_peaks columns per sample
    fripDF = pd.DataFrame(countsStackedDF.groupby(['sampleID','mapped_reads'], as_index=False)['reads_in_region'].sum())
    
    # Calculate FRiP in each sample
    fripDF['frip'] = (fripDF['reads_in_region'] / fripDF['mapped_reads']) * 100
    
    # Add number of regions for each (should be all the same number)
    fripDF['num_regions'] = countsStackedDF.groupby(['sampleID'])['apn'].count().unique()[0]

    # Print out FRiP TSV to user provided file
    fripDF.to_csv(args.outTSV, encoding='utf-8', index=False, sep="\t")
    print("FRiP table TSV saved to "+args.outTSV)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
