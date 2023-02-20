#!/usr/bin/env python3

import pandas as pd
import numpy as np
import scipy.stats as st
import math
import sqlite3
import argparse
import sys


## Use the following to get stacked with sampleID
## for file in cvrg_cnts_*; do 
##     NAME=$(basename $file _combined.csv | sed 's/cvrg_cnts_//')
##     awk -v name=$NAME '{if(NR==1){print "sampleID,"$0;}else{print name","$0}}' $file > temp_$NAME.csv
## done
## awk 'FNR==1 && NR!=1{next;}{print}' temp_* > combined_cnts.csv
## rm temp_*

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in a single stacked counts matrix (cat multiple together) to determine APN detection above background (DABG) and output TSV of flags")

    # Input data
    parser.add_argument("-i", "--input-stacked", dest="inCounts", required=True, help="Input CSV of counts in stacked format")
    parser.add_argument("-p", "--p-value", dest="pValue", required=False, type=float, help="P-value used in z-tests to flag DABG regions", default=0.05)
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temp SQL database file to reduce memory required", default=":memory:")
    parser.add_argument('-f', "--flag-reps", dest="flagReps", required=True, help="CSV file of flags for peaks called in each sample")

    # Output data
    parser.add_argument("-o", "--output", dest="outFlags", required=True, help="Output file name for TSV of DABG flags")
    parser.add_argument("-t", "--full-table", dest="outTable", default="", help="Output the full tall table of values used in DABG as a TSV to filename provided.")

    args = parser.parse_args()
    return args
def main():
    # Check that P value is valid
    if args.pValue > 1 or args.pValue < 0:
        sys.stderr.write("Invalid P value - must be between 0 and 1")
        sys.exit()

    # Get replicate flags for peaks
    repFlagsDF = pd.read_csv(args.flagReps)
    # Maybe could check format of flagReps file - check that all rows after coordinates are 0 or 1
    #if repFlagsDF.iloc[-1]
    repFlagsDF = repFlagsDF.rename(columns={'peak_id':'peakID'})
    for rep in range(4,len(list(repFlagsDF))) :
        repFlagsDF = repFlagsDF.rename(columns={repFlagsDF.columns[rep]:repFlagsDF.columns[rep].split('/')[-1].split("nomod147ext_")[1]})
    
    # Melt replicate flags into a stacked format so it can be merged later
    repFlagsTall= pd.melt(repFlagsDF.drop(['chrom','peak_start','peak_end'], axis=1), id_vars=['peakID'], var_name='sampleID', value_name='flag_peakCalled')
    
    # Connect to SQL database
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    
    # Import counts data and output sampleID names to stdout
    countsStackedDF = pd.read_csv(args.inCounts)
    countsStackedDF = countsStackedDF.rename(columns={'featureID':'peakID'})
    
    # Check that each sampleID has the same number of peaks
    count = len(countsStackedDF.groupby(['sampleID'])['apn'].count().unique())
    if count != 1 :
        sys.stderr.write("ERROR : Number of Peaks is not the same for all samples")
        sys.exit()
    numPeaks = countsStackedDF.groupby(['sampleID'])['apn'].count()[1]
    
    # Calculate Q1 and Std dev for each sampleID
    lowerQuartileDF=countsStackedDF.groupby(['sampleID'])['apn'].quantile(0.25).reset_index().rename(columns={'index':'sampleID', 'apn':'apnQ1'})
    stdDevDF = pd.DataFrame(countsStackedDF.groupby(['sampleID'])['apn'].std()).reset_index().rename(columns={'index':'sampleID', 'apn':'apnStdDev'})
    
    # Send counts dataset, called peak flag, lower quartiles, and standard deviations to SQL and merge
    # first merge: stdDev to lowerQuartile by sampleID
    # second merge: stdDev+lowerQuartile to countsTall by sampleID
    # third merge: stdDev+lowerQuartile+countsTall to flag_peakCalled by sampleID and peakID
    countsStackedDF.to_sql("countsTall", con, if_exists="replace")
    repFlagsTall.to_sql("flagPeak", con, if_exists="replace")
    lowerQuartileDF.to_sql("lowerQuartile", con, if_exists="replace")
    stdDevDF.to_sql("stdDev", con, if_exists="replace")
    cur.execute("CREATE TABLE statsKey AS SELECT in1.apnStdDev, in2.sampleID, in2.apnQ1 "
                "FROM stdDev in1 INNER JOIN lowerQuartile in2 "
                "ON in1.sampleID = in2.sampleID ;")
    cur.execute("CREATE TABLE countsKey AS SELECT in1.*, in2.peakID, in2.apn "
                "FROM statsKey in1 INNER JOIN countsTALL in2 "
                "ON in1.sampleID = in2.sampleID ;")
    cur.execute("CREATE TABLE flagKey AS SELECT in1.*, in2.flag_peakCalled "
                "FROM countsKey in1 INNER JOIN flagPeak in2 "
                "ON in1.sampleID = in2.sampleID and in1.peakID = in2.peakID ;")
    mergeDF = pd.read_sql("SELECT * FROM flagKey", con)
    
    # Calculate Z score and P value for each element
    sqrtN = round(math.sqrt(numPeaks))
    mergeDF['Zscore'] = (mergeDF['apn'] - mergeDF['apnQ1'])*sqrtN/mergeDF['apnStdDev']
    mergeDF['pValue'] = st.norm.sf(abs(mergeDF['Zscore']))
    
    # Determine whether P values are less than or equal to the user provided p value and set flag
    mergeDF['flag_DABG'] = np.where(mergeDF['pValue'] <= args.pValue, 1, 0)

    # Get name of sample type and split it into antibody, species, sex, and treatment
    name = "_".join(list(repFlagsDF)[4].split("_")[:-1])
    antibody = name.split("_")[0]
    species = name.split("_")[1]
    sex = name.split("_")[3]
    treatment = name.split("_")[4]
    
    # Output summary counts
    # sampleType,consensus_peaks,called_all_reps,not_called_rep1(#DABG),not_called_rep2(#DABG),not_called_rep3(#DABG)
    tempFlags = mergeDF.groupby('peakID')['flag_peakCalled'].sum().rename('sum_peakCalled')
    calledAll = tempFlags[tempFlags==3].count()
    notCalledRep1 = mergeDF.where((mergeDF['sampleID']==name+"_rep1") & (mergeDF['flag_peakCalled']==0)).count()[1]
    notCalledRep2 = mergeDF.where((mergeDF['sampleID']==name+"_rep2") & (mergeDF['flag_peakCalled']==0)).count()[1]
    notCalledRep3 = mergeDF.where((mergeDF['sampleID']==name+"_rep3") & (mergeDF['flag_peakCalled']==0)).count()[1]
    DABGnotCalledRep1 = mergeDF.where((mergeDF['sampleID']==name+"_rep1") & (mergeDF['flag_peakCalled']==0) & (mergeDF['flag_DABG']==1)).count()[1]
    DABGnotCalledRep2 = mergeDF.where((mergeDF['sampleID']==name+"_rep2") & (mergeDF['flag_peakCalled']==0) & (mergeDF['flag_DABG']==1)).count()[1]
    DABGnotCalledRep3 = mergeDF.where((mergeDF['sampleID']==name+"_rep3") & (mergeDF['flag_peakCalled']==0) & (mergeDF['flag_DABG']==1)).count()[1]
    DABGall = calledAll + mergeDF.where((mergeDF['flag_peakCalled']==0) & (mergeDF['flag_DABG']==1)).count()[1]

    print("{},{},{},{},{},{},{}({}),{}({}),{}({}),{}".format(antibody,species,sex,treatment,numPeaks,calledAll,notCalledRep1,DABGnotCalledRep1,notCalledRep2,DABGnotCalledRep2,notCalledRep3,DABGnotCalledRep3,DABGall))
    

    # merge called peak flag sum with full dataset
    fullMergeDF = mergeDF.merge(pd.DataFrame(tempFlags),how='left',on='peakID',validate="many_to_one")

    # Print out full table if requested by user
    if args.outTable != "" :
        fullMergeDF.to_csv(args.outTable, encoding='utf-8', index=False, sep="\t")

    # Print out flag file
    # The 'mean' is only used to aggregate the values (only one value so the average is that value)
    flagDF = mergeDF.groupby(['peakID','sampleID'])['flag_DABG'].aggregate('mean').unstack().reset_index()
    flagDF.to_csv(args.outFlags, encoding='utf-8', index=False, sep="\t")
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

