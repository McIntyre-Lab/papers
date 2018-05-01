#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-04
# NAME: add_exon_info_to_junctions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports the annotations for exons and adds them to the set of logical junctions
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
    parser = argparse.ArgumentParser(description="Import fragments and fusion annotation files")
    parser.add_argument("--input-junction-file", dest="inJuncInfo", action='store', required=True,
                        help="Input exon junction and exon-intron border junction CSV file")
    parser.add_argument("--input-exon-file", dest="inExonInfo", action='store', required=True,
                        help="Input exon info CSV file")
    parser.add_argument("--output-junction-info", dest="outCSV", action='store', required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args


def main():
    # Connect to database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Import concatenated set of exon-exon junctions and exon-intron border junctions
    juncInDF = pd.read_csv(args.inJuncInfo, sep=",", skiprows=1,
                           names=['event_id','event_type','chr','strand','num_transcripts','transcript_id',
                                  'donor_exon','acceptor_exon','donor_start','donor_stop','donor_size',
                                  'acceptor_start','acceptor_stop','acceptor_size','flag_junction_annotated',
                                  'flag_border_junction','flag_exonskip','num_skipped_exons','cat_skipped_exons'])
    # Import exon info
    exonInfoDF = pd.read_csv(args.inExonInfo, sep=",",
                             usecols=['gene_id', 'exon_id', 'exon_group', 'flag_alt_donor', 'flag_alt_acceptor'])

    juncInDF.to_sql("juncInfo", con, if_exists="replace")
    exonInfoDF.to_sql("exonInfo", con, if_exists="replace")

    # Merge donor and acceptor information. If donor is an intron, then this info is skipped
    cur.execute("CREATE TABLE donorInfo AS SELECT distinct gene_id AS donor_gene, exon_id AS donor_exon,"
                "exon_group AS donor_group, flag_alt_donor FROM exonInfo ORDER BY donor_exon ;")
    cur.execute("CREATE TABLE juncWdonor AS SELECT in1.*, in2.donor_gene, in2.donor_group, in2.flag_alt_donor "
                "FROM juncInfo in1 LEFT JOIN donorInfo in2 "
                "ON in1.donor_exon = in2.donor_exon ; ")

    cur.execute("CREATE TABLE acceptorInfo AS SELECT distinct gene_id AS acceptor_gene, exon_id AS acceptor_exon,"
                "exon_group AS acceptor_group, flag_alt_acceptor FROM exonInfo ORDER BY acceptor_exon ;")
    cur.execute("CREATE TABLE juncWacceptor AS SELECT in1.*, in2.acceptor_gene, in2.acceptor_group, in2.flag_alt_acceptor "
                "FROM juncWdonor in1 LEFT JOIN acceptorInfo in2 "
                "ON in1.acceptor_exon = in2.acceptor_exon ; ")
    juncWexonInfoDF = pd.read_sql("SELECT event_id, event_type, chr, strand, num_transcripts, transcript_id, donor_exon, "
                                  "acceptor_exon, donor_start, donor_stop, donor_size, acceptor_start, acceptor_stop, "
                                  "acceptor_size, flag_junction_annotated, flag_border_junction, flag_exonskip, "
                                  "num_skipped_exons, cat_skipped_exons, donor_gene, donor_group, flag_alt_donor, acceptor_gene, "
                                  "acceptor_group, flag_alt_acceptor FROM juncWacceptor", con )

    # Replace missing values with 0 (e.g. for alt donor/acceptor flags
    juncWexonInfoDF = juncWexonInfoDF.fillna({'num_transcripts':0,'donor_start':0,'donor_stop':0,'donor_group':0,
                                      'flag_alt_donor':0,'acceptor_start':0,'acceptor_stop':0,'acceptor_group':0,
                                      'flag_alt_acceptor':0,'num_skipped_exons':0,'transcript_id': 'Unannotated',
                                              'cat_skipped_exons':'(none)'})
    # Output annotation CSV
    with open(args.outCSV, 'w') as outCSVFile:
        juncWexonInfoDF.to_csv(outCSVFile, encoding='utf-8', index=False, sep=',')

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    main()
    logger.info('Script complete')
