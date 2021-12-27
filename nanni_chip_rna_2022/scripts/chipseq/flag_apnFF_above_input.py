#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sqlite3
import argparse


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in a single counts file to determine APN detection above input and output CSV of flags")

    # Input data
    parser.add_argument("-c", "--chipCounts", dest="chipCounts", required=True, help="CSV of ChIP (K4 or K27) counts from coverage counts step")
    parser.add_argument("-i", "--inputCounts", dest="inputCounts", required=True, help="CSV of corresponding Input counts from coverage counts step")
    parser.add_argument("-m", "--medianMapped", dest="medianMapped", required=True, help="Median mapped reads within the species")
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temp SQL database file to reduce memory required", default=":memory:")

    # Output data
    parser.add_argument("-o", "--output", dest="outFlags", required=True, help="Output file name for CSV of detection above input (DAI) flags")

    args = parser.parse_args()
    return args

def calc_ff(counts):
    # Calculate median "fudge factor"
    # ff = median mapped reads divided by mapped reads in the sample
    median = args.medianMapped
    mappedReads = counts['mapped_reads'].astype(int)[0]
    return float(median)/int(mappedReads)

def main():
        
    # Import ChIP and Input counts data
    chipDF = pd.read_csv(args.chipCounts)
    inputDF = pd.read_csv(args.inputCounts)

    # Calculate median mapped reads fudge factor
    chipFF = calc_ff(chipDF)
    inputFF = calc_ff(inputDF)

    # Multiply APN by ff
    chipDF['apnFF'] = chipDF['apn'] * chipFF
    inputDF['apnFF'] = inputDF['apn'] * inputFF

    # Connect to SQL database
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()

    # Merge input apnFF value into the chip dataframe
    chipDF.to_sql("chip", con, if_exists="replace")
    inputDF.rename(columns={"apnFF":"input_apnFF"}).to_sql("input", con, if_exists="replace")
    cur.execute("CREATE TABLE merged AS SELECT in1.featureID, in1.apnFF, in2.input_apnFF "
                "FROM chip in1 LEFT JOIN input in2 "
                "ON in1.featureID = in2.featureID ;")
    mergeDF = pd.read_sql("SELECT * FROM merged", con)

    # Flag detection above input (DAI)
    mergeDF['flag_dai'] = np.where(mergeDF['apnFF'] > mergeDF['input_apnFF'], 1, 0)

    # Output flags
    mergeDF.to_csv(args.outFlags, index=False)

   
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

