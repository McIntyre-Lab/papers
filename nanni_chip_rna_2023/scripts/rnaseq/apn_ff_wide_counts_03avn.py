#!/usr/bin/env python

import pandas as pd
import sqlite3
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Multiply sampleID specific FF to wide counts")

    # Input data
    parser.add_argument("-f", "--inFF", dest="inFF", required=True, help="Input FF TSV file with sampleID column")
    parser.add_argument("-n", "--inName", dest="inName", required=True, help="FF column name")
    parser.add_argument("-c", "--inCount", dest="inCount", required=True, help="Input wide format APN coverage count file")

    # Output data
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temporary SQL database file to reduce memory requiredfor SQL merge", default=":memory:")
    parser.add_argument("-w", "--outWide", dest="outWide", required=False, help="Output CSV file name for wide counts multiplied by FF")
    parser.add_argument("-s", "--outStack", dest="outStack", required=False, help="Output CSV file name for stacked counts multiplied by FF")

    args = parser.parse_args()
    return args

def sexID(s):
    if s.count("_") == 4:
        return s.split("_")[2]
    else:
        return s

def main():

    # Get FF file and column name
    ffDF = pd.read_csv(args.inFF, sep="\t")
    name = args.inName

    # Drop any columns that are flag_bad_sample
    ffBAD = ffDF[ffDF['flag_bad_sample']==1]['sampleID'].tolist()
    ffDF = ffDF[ffDF['flag_bad_sample']==0]

    # Get coverage count file
    countDF = pd.read_csv(args.inCount).set_index('featureID')
    
    # Drop any bad samples
    cols = [c for c in countDF.columns if c not in ffBAD]
    countDF = countDF[cols]

    # Stack coverage counts (flags dropped)
    stackDF = countDF[[c for c in countDF.columns if "flag" not in c]].stack().reset_index(
              ).rename(columns={'level_1':'sampleID',0:'apn'})

    # Merge FF and counts by sampleID
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()

    ffDF.to_sql("FF", con, if_exists="replace")
    stackDF.to_sql("count", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2."+name+" "
                "FROM count in1 LEFT JOIN FF in2 "
                "ON in1.sampleID = in2.sampleID ")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])
    
    # Multiply apn by FF
    mergeDF['apn_ff'] = mergeDF['apn'] * mergeDF[name]

    # Output stacked if requested (flags will not be included)
    if args.outStack:
        mergeDF.to_csv(args.outStack,index=False)
    
    # Reformat to wide
    wideDF = mergeDF.pivot(index='featureID',columns='sampleID',values='apn_ff')
    wideDF = wideDF.reindex(sorted(wideDF.columns,key=sexID),axis=1)
    wideDF.columns = wideDF.columns + "_apn_uq_ff"
    
    # Merge normalized APN with raw APN
    countDF.columns = countDF.columns + "_apn"
    countDF = countDF.rename(columns={'featureID_apn':'featureID'})
    countDF.reset_index().to_sql("raw", con, if_exists="replace")
    wideDF.reset_index().to_sql("norm", con, if_exists="replace")
    cur.execute("CREATE TABLE full AS SELECT in1.*, in2.* "
                "FROM raw in1 LEFT JOIN norm in2 "
                "ON in1.featureID = in2.featureID ")
    fullDF = pd.read_sql("SELECT * FROM full", con).drop(columns=['index','index:1','featureID:1'])

    # Output wide of all counts (raw and normalized) if requested
    if args.outWide:
        fullDF.to_csv(args.outWide, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

