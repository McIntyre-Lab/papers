#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-09
# NAME: append_border_junctions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: Append border junctions (sequences that span exon-intron borders) to annotated junction catalog and
#              add missing information
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
    parser.add_argument("--input-donor-border-junctions", dest="inDonorBorders", action='store',
                        required=True, help="Input donor exon-intron border junctions")
    parser.add_argument("--input-acceptor-border-junctions", dest="inAcceptorBorders", action='store',
                        required=True, help="Input acceptor exon-intron border junctions")
    parser.add_argument("--junction-size", dest="juncSize", action='store', type=int, required=True, help="Junction size")
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
    logger.info('Importing exon-intron border junctions')
    donorBordersDF = pd.read_csv(args.inDonorBorders, sep="\t",
                          names=['chr','totalstart','totalstop','event_id','score','strand','color',
                                 'blocks','block_sizes','block_starts'])
    acceptorBordersDF = pd.read_csv(args.inAcceptorBorders, sep="\t",
                          names=['chr','totalstart','totalstop','event_id','score','strand','color',
                                 'blocks','block_sizes','block_starts'])
    bordersDF = pd.concat([donorBordersDF,acceptorBordersDF])
    # Put dataframes into SQL
    juncDF.to_sql("juncInfo", con, if_exists="replace")
    bordersDF.to_sql("borderJunc", con, if_exists="replace")
    # Format border junction info, split into donor/acceptor pieces
    cur.execute("CREATE TABLE borderJunc2 AS SELECT chr, event_id, strand, substr(event_id, 1, exPos-1) AS donor_exon, "
                "substr(event_id, exPos+1) AS acceptor_exon, totalstart AS donor_start, totalstop AS acceptor_stop, block_sizes "
                "FROM  (SELECT *, instr(event_id,'|') AS exPos FROM borderJunc ) "
                "ORDER BY event_id; ")
    # Add some varibles: donor_stop, acceptor_start, donor_size, acceptor_size
    #                    transcript_id ('Unannotated'), num_transcripts, flag_junction_annotated,flag_exonskip,
    #                    num_skipped_exons, cat_skipped_exons, event_type ("exon_intron_border")
    cur.execute("ALTER TABLE borderJunc2 ADD donor_size INT ;")
    cur.execute("UPDATE borderJunc2 SET donor_size = (CASE WHEN strand='+' AND donor_exon='intron' THEN (? - 1) "
                "WHEN strand='+' AND acceptor_exon='intron' THEN (block_sizes - ?) "
                "WHEN strand='-' AND donor_exon='intron' THEN ? "
                "WHEN strand='-' AND acceptor_exon='intron' THEN (block_sizes - ? - 1) end); ",
                (args.juncSize,args.juncSize,args.juncSize,args.juncSize))
    cur.execute("ALTER TABLE borderJunc2 ADD acceptor_size INT ;")
    cur.execute("UPDATE borderJunc2 SET acceptor_size = (CASE WHEN strand='+' AND donor_exon='intron' THEN (block_sizes - ?) "
                "WHEN strand='+' AND acceptor_exon='intron' THEN (? - 1) "
                "WHEN strand='-' AND donor_exon='intron' THEN (block_sizes - ? - 1) "
                "WHEN strand='-' AND acceptor_exon='intron' THEN ? end); ",
                (args.juncSize,args.juncSize,args.juncSize,args.juncSize))
    cur.execute("ALTER TABLE borderJunc2 ADD donor_stop INT ;")
    cur.execute("UPDATE borderJunc2 SET donor_stop = donor_start + donor_size ;")
    cur.execute("ALTER TABLE borderJunc2 ADD acceptor_start INT ;")
    cur.execute("UPDATE borderJunc2 SET acceptor_start = acceptor_stop - acceptor_size ;")
    cur.execute("ALTER TABLE borderJunc2 ADD transcript_id TEXT ;")
    cur.execute("UPDATE borderJunc2 SET transcript_id = 'Unannotated' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD num_transcripts INT ;")
    cur.execute("UPDATE borderJunc2 SET num_transcripts = '0' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD flag_junction_annotated INT ;")
    cur.execute("UPDATE borderJunc2 SET flag_junction_annotated = '0' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD flag_exonskip INT ;")
    cur.execute("UPDATE borderJunc2 SET flag_exonskip = '0' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD num_skipped_exons INT ;")
    cur.execute("UPDATE borderJunc2 SET num_skipped_exons = '0' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD cat_skipped_exons TEXT ;")
    cur.execute("UPDATE borderJunc2 SET cat_skipped_exons = '' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD event_type TEXT ;")
    cur.execute("UPDATE borderJunc2 SET event_type = 'exon_intron_border' ;")
    cur.execute("ALTER TABLE borderJunc2 ADD flag_border_junction INT ;")
    cur.execute("UPDATE borderJunc2 SET flag_border_junction = '1' ;")
    # Append exon-exon junctions

    cur.execute("ALTER TABLE juncInfo ADD flag_border_junction INT ;")
    cur.execute("UPDATE juncInfo SET flag_border_junction = '0' ;")

    cur.execute("CREATE TABLE juncAll AS SELECT * FROM borderJunc2 ; ")
    cur.execute("INSERT INTO  juncAll(chr, event_id, strand, donor_exon, donor_size, donor_start, "
                "donor_stop, acceptor_exon, acceptor_size, acceptor_start, acceptor_stop, event_type, num_transcripts, "
                "transcript_id, flag_junction_annotated, flag_exonskip, num_skipped_exons, cat_skipped_exons, flag_border_junction) "
                "SELECT chr, event_id, strand, donor_exon, donor_size, donor_start, "
                "donor_stop, acceptor_exon, acceptor_size, acceptor_start, acceptor_stop, event_type, num_transcripts, "
                "transcript_id, flag_junction_annotated, flag_exonskip, num_skipped_exons, cat_skipped_exons, flag_border_junction "
                "FROM juncInfo ; ")
    juncAllDF=pd.read_sql("SELECT * FROM juncAll ORDER BY event_id;", con)
    juncAllDF = juncAllDF[['event_id','event_type','chr','strand','num_transcripts','transcript_id','donor_exon',
                  'acceptor_exon','donor_start','donor_stop','donor_size','acceptor_start','acceptor_stop',
                  'acceptor_size','flag_junction_annotated','flag_border_junction','flag_exonskip',
                  'num_skipped_exons','cat_skipped_exons']]
    logger.info("Exporting exon-exon junctions and exon-intron border junctions")
    with open(args.outCSV, 'w') as outFile:
        juncAllDF.to_csv(outFile, encoding='utf-8', index=False)

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