#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sqlite3
import argparse


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in flags of detection above input in all replicates to determine presence (at least 2 of 3 replicates DAI) and output CSV of flags")

    # Input data
    parser.add_argument("-n","--name", dest="inName", required=True, help="Name of sample type")
    parser.add_argument("-1", "--rep1", dest="inREP1", required=True, help="Replicate 1 CSV of detection above input (DAI) flags")
    parser.add_argument("-2", "--rep2", dest="inREP2", required=True, help="Replicate 2 CSV of detection above input (DAI) flags")
    parser.add_argument("-3", "--rep3", dest="inREP3", required=True, help="Replicate 3 CSV of detection above input (DAI) flags")
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temp SQL database file to reduce memory required", default=":memory:")

    # Output data
    parser.add_argument("-o", "--output", dest="outFlags", required=True, help="Output file name for CSV of presence of flags")

    args = parser.parse_args()
    return args

def calc_ff(counts):
    # Calculate median "fudge factor"
    # ff = median mapped reads divided by mapped reads in the sample
    median = args.medianMapped
    mappedReads = counts['mapped_reads'][0]
    return float(median)/int(mappedReads)

def main():
    
    # Get name of sample type
    name=args.inName
    
    # Import DAI flag files
    Rep1DF = pd.read_csv(args.inREP1).rename(columns={'apnFF':'rep1_apnFF', 'input_apnFF':'rep1_input_apnFF', 'flag_dai':'flag_dai_rep1'})
    Rep2DF = pd.read_csv(args.inREP2).rename(columns={'apnFF':'rep2_apnFF', 'input_apnFF':'rep2_input_apnFF', 'flag_dai':'flag_dai_rep2'})
    Rep3DF = pd.read_csv(args.inREP3).rename(columns={'apnFF':'rep3_apnFF', 'input_apnFF':'rep3_input_apnFF', 'flag_dai':'flag_dai_rep3'})


    # Connect to SQL database
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()

    # Merge input apnFF value into the chip dataframe
    Rep1DF.to_sql("rep1", con, if_exists="replace")
    Rep2DF.to_sql("rep2", con, if_exists="replace")
    cur.execute("CREATE TABLE first AS SELECT in1.featureID, in1.rep1_apnFF, in1.rep1_input_apnFF, in1.flag_dai_rep1, in2.rep2_apnFF, in2.rep2_input_apnFF, in2.flag_dai_rep2 "
                "FROM rep1 in1 LEFT JOIN rep2 in2 "
                "ON in1.featureID = in2.featureID ;")
    Rep3DF.to_sql("rep3", con, if_exists="replace")
    cur.execute("CREATE TABLE second AS SELECT in1.*, in2.rep3_apnFF, in2.rep3_input_apnFF, in2.flag_dai_rep3 "
                "FROM first in1 LEFT JOIN rep3 in2 "
                "ON in1.featureID = in2.featureID ;")  
    mergeDF = pd.read_sql("SELECT * FROM second", con)

    # Flag detection above input (DAI)
    mergeDF['flag_presence'] = np.where(mergeDF['flag_dai_rep1']+mergeDF['flag_dai_rep2']+mergeDF['flag_dai_rep3'] >= 2, 1, 0)

    # Output flags
    mergeDF.to_csv(args.outFlags, index=False)
    
    # Calculate summary
    totalFeat = mergeDF.count()[1]
    allDAI = mergeDF.where(mergeDF['flag_dai_rep1']+mergeDF['flag_dai_rep2']+mergeDF['flag_dai_rep3'] == 3).count()[1]
    notDAIrep1 = mergeDF.where((mergeDF['flag_dai_rep1']==0) & (mergeDF['flag_dai_rep2']+mergeDF['flag_dai_rep3']==2)).count()[1]
    notDAIrep2 = mergeDF.where((mergeDF['flag_dai_rep2']==0) & (mergeDF['flag_dai_rep1']+mergeDF['flag_dai_rep3']==2)).count()[1]
    notDAIrep3 = mergeDF.where((mergeDF['flag_dai_rep3']==0) & (mergeDF['flag_dai_rep1']+mergeDF['flag_dai_rep2']==2)).count()[1]

    print("{},{},{},{},{},{}".format(name,totalFeat,allDAI,notDAIrep1,notDAIrep2,notDAIrep3))


   
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

