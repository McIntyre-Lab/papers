#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2017-12-19
# NAME: import_and_format_junctions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports the junctions BED file created by extractJunctions.py and reformats it for adding
# in donor and acceptor exon information
#
# REQUIRED PACKAGES: argparse  (tested with v1.1)
#                    pandas    (tested with v0.19.2)
#                    logging   (tested with v0.5.1.2)
#
#######################################################################################################################

# Packages

import argparse
import logging
import sqlite3
import pandas as pd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import junctions BED file")
    parser.add_argument("--bed", dest="bedInput", action='store', required=True, help="Input BED file")
    parser.add_argument("--output", dest="outCSV", action='store', required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args

def main():
    # Import exon, fragment and junction annotations
    logger.info("Importing BED file of all possible, logical exon-exon junctions")
    juncDF = pd.read_csv(args.bedInput, sep="\t", names=['chr','totalStart','totalStop','event_id','score','strand',
                                                         'totalStart2','totalStop2','color','blocks','block_sizes',
                                                         'block_starts'])

    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Convert dataframe to SQL table
    cur.execute("CREATE TABLE juncInfo (chr TEXT, totalStart INT, totalStop INT,event_id TEXT, score INT, strand TEXT, "
                "totalStart2 INT, totalStop2 INT, color TEXT, blocks  INT, block_sizes TEXT, block_starts TEXT) ;")
    juncDF.to_sql("juncInfo", con, if_exists="replace")
    # Get donor and acceptor exon info by splitting event_id, block_sizes, and calculate start/stop positions
    cur.execute("CREATE TABLE juncInfo2 AS SELECT chr, strand, event_id, substr(event_id, 1, exPos-1) AS donor_exon, "
                "substr(block_sizes, 1, exSize-1) AS donor_size, totalStart AS donor_start, "
                "(totalStart + substr(block_sizes, 1, exSize-1)) AS donor_stop, "
                "substr(event_id, exPos+1 ) AS acceptor_exon, substr(block_sizes, exSize+1) AS acceptor_size, "
                "(totalStop - substr(block_sizes, exSize+1)) AS acceptor_start, totalStop AS acceptor_stop "
                "FROM (SELECT *, instr(event_id,'|') AS exPos, instr(block_sizes,',') AS exSize FROM juncInfo ) "
                "ORDER BY event_id ;")
    cur.execute("ALTER TABLE juncInfo2 ADD event_type text ;")
    cur.execute("UPDATE juncInfo2 SET event_type = 'exon_junction' ;")

    juncDF2 = pd.read_sql("SELECT * from juncInfo2", con)

    # Export CSV of formatted junctions
    logger.info("Export formatted exon junctions")
    with open(args.outCSV, 'w') as outFile:
        juncDF2.to_csv(outFile, encoding='utf-8', index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    #print("Starting script...")
    main()
    logger.info('Script complete')
    #print("Script done!")