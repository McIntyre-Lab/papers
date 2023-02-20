#!/usr/bin/env python

import pandas as pd
import sqlite3
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge on/off counts and flags")

    # Input data
    parser.add_argument("-1", "--in1", dest="in1", required=True, help="First input fragment on/off flag CSV file")
    parser.add_argument("-2", "--in2", dest="in2", required=True, help="Second input fragment on/off flag CSV file (only flags will be added to first input)")
    parser.add_argument("-c", "--counts", dest="inCounts", required=True, help="Merged counts file")

    # Output data
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temporary SQL database file to reduce memory requiredfor SQL merge", default=":memory:")
    parser.add_argument("-o", "--outFile", dest="outFile", required=False, help="Output CSV file name for merged flags")

    args = parser.parse_args()
    return args

def main():

    # Get input flag files and counts
    firstDF = pd.read_csv(args.in1)
    secondDF = pd.read_csv(args.in2)
    countDF = pd.read_csv(args.inCounts)

    # Merge flags by featureID
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()

    firstDF.to_sql("first", con, if_exists="replace")
    secondDF.to_sql("second", con, if_exists="replace")
    countDF.to_sql("counts", con, if_exists="replace")
    cur.execute("CREATE TABLE flags AS SELECT in1.*, in2.* "
                "FROM first in1 LEFT JOIN second in2 "
                "ON in1.featureID = in2.featureID ")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.* "
                "FROM flags in1 LEFT JOIN counts in2 "
                "ON in1.featureID = in2.featureID ")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index','index:1','index:2','featureID:1','featureID:2'])
    
    mergeDF.to_csv(args.outFile,index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

