#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-04-13
# NAME: import_counts_and_convert_to_wide.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports all coverage count CSVs (generated from rpkm_calculate.py) in a folder and converts
# them from their native tall/long format into a single wide format CSV, where columns are samples and rows are
# events. Note that this scripts assumes the counts CSVs are in the same format as the output generated from
# rpkm_calculate.py.
# #
# REQUIRED PACKAGES: pandas    (tested with v0.19.2)
#                    argparse  (tested with v1.1)
#                    logging   (tested with v0.5.1.2)
#                    sqlite3
#                    glob
#
#######################################################################################################################

# Import required packages
import pandas as pd
import logging
import glob
import os
import sqlite3
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes a directory of coverage counts for events (generated"
                                                 "with rpkm_calculate.py and outputs a single CSV in wide format")
    parser.add_argument("--input-directory", dest="inPath", required=True, help="Path to coverage counts files")
    parser.add_argument("--output-wide", dest="outCSV", required=True, help="Output CSV of counts in wide format")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Iterate over each file in the input path and append to countsTall table
    counter=1
    for filename in glob.glob(os.path.join(args.inPath, '*.csv')):
        # If first file in list, then this dataset become the "base" dataset into which we will merge all
        # subsequent counts
        countsDF = pd.read_csv(filename, sep=",", skiprows=1,
                               names=['sample_id','event_id','mapped_reads','read_length','region_length',
                                      'region_depth','reads_in_region','apn','rpkm','mean','std','cv'],
                               usecols=['sample_id','event_id','apn'])
        sampleID=countsDF['sample_id'][0]
        countsDF.to_sql("countsTemp", con, if_exists="replace")
        if counter == 1:
            cur.execute("CREATE TABLE countsWide AS SELECT event_id, apn AS "+sampleID+" FROM countsTemp ORDER BY event_id;")
            counter=counter+1
        else :
            cur.execute("CREATE TABLE countsWide2 AS SELECT in1.*, in2.apn AS "+sampleID+" "
                            "FROM countsWide in1 INNER JOIN countsTemp in2 "
                            "ON in1.event_id = in2.event_id  ORDER BY event_id ;" )
            cur.execute("DROP TABLE countsWide;")
            cur.execute("ALTER TABLE countsWide2 RENAME TO countsWide ;")
    countsWideDF = pd.read_sql("SELECT * FROM countsWide;", con)

    # Write output index
    with open(args.outCSV, 'w') as outCounts:
        countsWideDF.to_csv(outCounts, encoding='utf-8', index=False, sep="\t")


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()

    # Setting up logger
    logger = logging.getLogger()
    logger.info('Starting script')

    # Calling main script
    main()
    logger.info('Script complete: index created!')


