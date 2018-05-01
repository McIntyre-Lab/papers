#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-05
# NAME: import_and_format_junctions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script flags which of the set of all possible, logical junctions derived from GFF are also
# annotated to known transcripts ("transcript-annotated junctions").
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
    parser = argparse.ArgumentParser(description="Flag transcript-annotated junctions")
    parser.add_argument("--input-junctions", dest="inJunc", action='store', required=True,
                        help="Input formatted logical junctions")
    parser.add_argument("--input-annotated-junctions", dest="inAnnotJunc", action='store',
                        required=True, help="Input transcript-annotated junctions")
    parser.add_argument("--output", dest="outCSV", action='store', required=True,help="Output CSV file name")
    args = parser.parse_args()
    return args

def main():
    # Connect to database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Import files
    logger.info('Importing formatted junctions')
    juncDF = pd.read_csv(args.inJunc, sep=",")
    logger.info('Importing transcript-annotated junctions')
    annotDF = pd.read_csv(args.inAnnotJunc, sep=",", skiprows=1,
                          names=['event_id','junc_coord','transcript_id','gene_id'],
                          usecols=['event_id','transcript_id'])
    # Convert to SQL tables
    cur.execute("CREATE TABLE juncInfo (chr TEXT, event_id TEXT, strand TEXT, acceptor_exon TEXT, acceptor_size INT, "
                "acceptor_start INT, acceptor_stop INT, donor_exon TEXT, donor_size INT, donor_start INT,"
                "donor_stop INT, event_type TEXT) ;")
    cur.execute("CREATE TABLE annotInfo (event_id TEXT, transcript_id TEXT) ;")
    juncDF.to_sql("juncInfo", con, if_exists="replace")
    annotDF.to_sql("annotInfo", con, if_exists="replace")

    # Count and concatenate transcripts by junction. Use "event_id" to merge on
    cur.execute("CREATE TABLE annotList AS SELECT distinct event_id, transcript_id "
                "FROM annotInfo ORDER BY event_id, transcript_id ;")
    cur.execute("CREATE TABLE annotList2 AS SELECT distinct event_id, count(distinct transcript_id) AS num_transcripts, "
                "replace(group_concat(distinct transcript_id),',', '|') AS transcript_id "
                "FROM annotList "
                "GROUP BY event_id ;")
    ## Flag annotated/unannotated junctions
    cur.execute("SELECT * FROM annotList2 ORDER BY event_id; ")
    cur.execute("SELECT * FROM juncInfo ORDER BY event_id; ")
    cur.execute("CREATE TABLE junc2annot AS SELECT in1.*, in2.transcript_id, in2.num_transcripts "
                "FROM juncInfo in1 LEFT JOIN annotList2 in2 "
                "ON in1.event_id = in2.event_id ;")
    #cur.execute("SELECT * FROM junc2annot")
    #allGene2Frag = cur.fetchall()
    #print(allGene2Frag)
    cur.execute("UPDATE junc2annot SET num_transcripts=0 WHERE coalesce(transcript_id, '') = '' ;")
    cur.execute("UPDATE junc2annot SET transcript_id='Unannotated' WHERE coalesce(transcript_id, '') = '' ;")
    cur.execute("ALTER TABLE junc2annot ADD flag_junction_annotated int ;")
    cur.execute("UPDATE junc2annot "
                "SET flag_junction_annotated = (case when num_transcripts=0 then '0' "
                "when num_transcripts > 0 then '1' end);")

    # Check output
    junc2annotDF=pd.read_sql("SELECT chr, event_id, strand, acceptor_exon, acceptor_size, acceptor_start, acceptor_stop, "
                             "donor_exon, donor_size, donor_start, donor_stop, event_type, num_transcripts, "
                             "transcript_id, flag_junction_annotated FROM junc2annot ORDER BY event_id ;", con)

    # Export
    logger.info("Exporting junctions with transcript info")
    with open(args.outCSV, 'w') as outFile:
        junc2annotDF.to_csv(outFile, encoding='utf-8', index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    # print("Starting script...")
    main()
    logger.info('Script complete')
    # print("Script done!")