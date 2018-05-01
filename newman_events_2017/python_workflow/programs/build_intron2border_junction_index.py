#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2017-12-15
# NAME: build_Event2Transcript_index.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script creates an intron-to-border junction index file used by Event Analysis to report
# the read coverage of introns, their associated border junctions and flanking exonic regions (fusions), to aid
# the user in deciding whether there is evidence on intron retention, alternative/novel splice usage, etc.
# It takes the annotation CSVs for junctions, exonic regions and introns to assemble a complete intron/border index,
# where each border junction and intron are assigned to a single intron event, flanked by its neighboring
# exonic regions. Where the exonic regions of intron events can be assigned to multiple genes, then the output of this
# intron event is suppressed, to avoid instances of overlapping intron events.
#
# REQUIRED PACKAGES: pandas    (tested with v0.19.2)
#                    argparse  (tested with v1.1)
#                    logging   (tested with v0.5.1.2)
#
#######################################################################################################################

# Import required packages
import pandas as pd
import logging
import argparse
import sqlite3

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Generate an intron-to-border-junction index file for"
                                                 "interpreting read coverage of intronic regions")
    parser.add_argument('--intron-annotation-file', dest="inIntrons", required=True, help="Input intron annotation CSV")
    parser.add_argument("--junction-annotation-file", dest="inJunctions", required=True,
                        help="Input junction annotation CSV")
    parser.add_argument("--output-intron-index", dest="outCSV", required=True,
                        help="Output event index CSV")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    # Import intron and junction annotations
    logger.info("Importing intron and junction annotations")
    intronDF = pd.read_csv(args.inIntrons, usecols=('intron_id','chr','intron_start','intron_stop','gene_id',
                                                    'exonic_region_id_5prime','exonic_region_id_3prime'))
    juncDF = pd.read_csv(args.inJunctions, usecols=('junction_id','chr','donor_stop','acceptor_start','transcript_id',
                                                       'gene_id','flag_border_junction'))
    # Convert to SQL tables
    intronDF.to_sql("intronInfo", con, if_exists="replace")
    juncDF.to_sql("juncInfo", con, if_exists="replace")

    # So border junctions and introns can be merged, donor_stop and acceptor start need to renamed to intron_start
    # and intron_stop respectively. When the "donor exon" is an intron, donor_stop = intron_stop
    # When the "acceptor exon" is an intron, acceptor_start = intron_start
    # I am going to first map 5' border junctions to the 5' end  of introns, then 3'
    # border junctions for the 3' end of the introns.
    # First, I want to expand concatenated gene IDs. Junctions with multiple gene ID shouldn't be retained in the
    # final output, but iterate over these for completeness
    cur.execute("""Select junction_id, chr , donor_stop , acceptor_start , gene_id from juncInfo WHERE flag_border_junction = 1""")
    allBorders = cur.fetchall()
    cur.execute("""CREATE TABLE IF NOT EXISTS borderInfo
                        (junction_id TEXT, chr TEXT, donor_stop INT, acceptor_start INT, gene_id TEXT)""")
    for border in allBorders:
        genes = border[4].split("|")
        for gn in genes:
            cur.execute("INSERT INTO borderInfo VALUES(:junction_id, :chr, :donor_stop, :acceptor_start, :gene_id)",
                           {"junction_id": border[0], "chr": border[1], "donor_stop": border[2],
                            "acceptor_start": border[3], "gene_id":gn})
    # Merge INNER with intron table on chromosome, gene, and acceptor_start (as intron_start)
    cur.execute("CREATE TABLE intronWstart AS SELECT in1.intron_id, in1.chr, in1.intron_start, in1.intron_stop, "
                "in1.gene_id, in1.exonic_region_id_5prime, in2.junction_id AS border_junction_id_5prime "
                "FROM intronInfo in1 INNER JOIN borderInfo in2 "
                "ON in1.chr = in2.chr AND in1.gene_id = in2.gene_id AND in1.intron_start = in2.acceptor_start ;")
    # Merge INNER with intron table on chromosome, gene, and donor_stop (as intron_stop)
    cur.execute("CREATE TABLE intronWstop AS SELECT in1.intron_id, in1.chr, in1.gene_id, "
                "in1.exonic_region_id_3prime, in2.junction_id AS border_junction_id_3prime "
                "FROM intronInfo in1 INNER JOIN borderInfo in2 "
                "ON in1.chr = in2.chr AND in1.gene_id = in2.gene_id AND in1.intron_stop = in2.donor_stop ;")
    cur.execute("CREATE TABLE intronBorderIndex AS SELECT in1.*, in2.exonic_region_id_3prime,"
                "in2.border_junction_id_3prime FROM intronWstart in1 "
                "INNER JOIN intronWstop in2 ON in1.gene_id = in2.gene_id AND in1.intron_id = in2.intron_id ;")

    intronBorderIndexDF = pd.read_sql("SELECT * FROM intronBorderIndex ORDER BY chr, intron_start, intron_stop ;", con)

    # Write output index
    with open(args.outCSV, 'w') as outIndex:
        intronBorderIndexDF.to_csv(outIndex, encoding='utf-8', index=False)

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