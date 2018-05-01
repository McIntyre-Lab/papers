#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-09
# NAME: import_exon_skipping_annotations.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: Import exon-skipping annotations for each logical junction and add to junction annotations
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
    parser.add_argument("--input-exonskip-annot", dest="inExonSkipAnnot", action='store',
                        required=True, help="Input exon-skipping annotations")
    parser.add_argument("--output", dest="outCSV", action='store', required=True,help="Output CSV file name")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Import data
    logger.info('Importing formatted junctions')
    juncDF = pd.read_csv(args.inJunc, sep=",")
    logger.info('Importing exon-skipping annotations')
    exonskipDF = pd.read_csv(args.inExonSkipAnnot, sep=",", skiprows=1,
                          names=['event_id','flag_exonskip','num_skipped_exons','cat_skipped_exons'])
    juncDF.to_sql("juncAnnot", con, if_exists="replace")
    exonskipDF.to_sql("exonSkip", con, if_exists="replace")
    cur.execute("CREATE TABLE juncExonSkip AS SELECT in1.*, in2.flag_exonskip, in2.num_skipped_exons, in2.cat_skipped_exons "
                "FROM juncAnnot in1 INNER JOIN exonSkip in2 "
                "ON in1.event_id = in2.event_id ;")
    cur.execute("UPDATE juncExonSkip SET cat_skipped_exons='' WHERE coalesce(cat_skipped_exons, '') = '' ;")
    juncExonSkipDF = pd.read_sql("SELECT chr, event_id, strand, donor_exon, donor_size, donor_start, donor_stop, "
                                 "acceptor_exon, acceptor_size, acceptor_start, acceptor_stop, event_type, num_transcripts, "
                                 "transcript_id, flag_junction_annotated, flag_exonskip, num_skipped_exons, "
                                 "cat_skipped_exons from juncExonSkip", con)
    # Export CSV
    logger.info("Exporting junctions with exon-skipping annotations info")
    with open(args.outCSV, 'w') as outFile:
        juncExonSkipDF.to_csv(outFile, encoding='utf-8', index=False)

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
