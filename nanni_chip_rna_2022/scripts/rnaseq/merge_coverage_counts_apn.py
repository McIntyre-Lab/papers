#!/usr/bin/env python

### Merge all coverage count files

import pandas as pd
import os
import logging
import argparse
import sqlite3

def getOptions():
    parser = argparse.ArgumentParser(description="Merge all coverage count files together (using apn values)")

    #Inputs
    parser.add_argument("--directory", "-d", dest="directory",action="store", required=True, help="Input directory of coverage files")
    parser.add_argument("--output", "-o", dest="output", action="store",required=True, help="Output file for merged coverage counts")

    args = parser.parse_args()
    return args

def main():

    # Change into input directory
    workingDIR = os.getcwd()
    os.chdir(args.directory)

    # Import coverage count files while merging
    num = 0
    mergeDF = pd.DataFrame()
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    for file in files :
        if num == 0:
            singleDF = pd.read_csv(file, sep=",", dtype='unicode')
            sampleID = "_".join(file.split(".")[0].split("_")[2:-1])
            mergeDF = singleDF.rename(columns={'apn':sampleID})[['featureID', sampleID]]
            num = num + 1
        else:
            singleDF = pd.read_csv(file, sep=",", dtype='unicode')
            sampleID = "_".join(file.split(".")[0].split("_")[2:-1])
            tempSingle = singleDF.rename(columns={'apn':sampleID})[['featureID',sampleID]]

            # Open new sqlite connection
            con = sqlite3.connect(":memory:")
            cur = con.cursor()
            tempSingle.sort_values('featureID').to_sql("newSample", con, if_exists="replace")
            mergeDF.sort_values('featureID').to_sql("allSample", con, if_exists="replace")

            # Inner join with other apn values
            cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.* "
                        "FROM allSample in1 INNER JOIN newSample in2 "
                        "ON in1.featureID = in2.featureID ")
            tempMerge = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index','index:1','featureID:1'])
            con.close()

            # Check lengths of dataframes for proper merge
            if len(tempSingle) != len(tempMerge) or len(mergeDF) != len(tempMerge):
                print("ERROR: IMPROPER MERGE WITH "+sampleID)
                exit()
            mergeDF = tempMerge
            num = num + 1

    # Output final merged file
    os.chdir(workingDIR)
    mergeDF.to_csv(args.output,index=False,header=True)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
