#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-04
# NAME: import_and_format_fragments.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports the annotations for exon fusions (exonic regions) extracted from GFF
# and formats them in preparation for creating event-to-transcript index
#
# REQUIRED PACKAGES: argparse  (tested with v1.1)
#                    pandas    (tested with v0.19.2)
#                    logging   (tested with v0.5.1.2)
#                    sqlite3   (tested with v2.6.0)
#
#######################################################################################################################

# Packages
import argparse
import logging
import sqlite3
import datetime
import pandas as pd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import fragments and fusion annotation files")
    parser.add_argument("--input-fusion-file", dest="inFusInfo", action='store', required=True,
                        help="Input fusion info TSV file")
    parser.add_argument("--input-fusion-bed", dest="inFusBED", action='store', required=True,
                        help="Input fusion BED file")
    parser.add_argument("--input-exon-file", dest="inExonInfo", action='store', required=True,
                        help="Input exon info CSV file")
    parser.add_argument("--outCSV", dest="outCSV", action='store', required=True, help="Output CSV file name -- this is the"
                                                                                       "complete annotation for exonic regions")
    parser.add_argument("--outBED", dest="outBED", action='store', required=True, help="Output BED file name -- this BED file"
                                                                                       "will be used for coverage counts")

    args = parser.parse_args()
    return args

def main():
    # Set up temporary SQL database in memory
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    # Import exon annotations and format
    logger.info("Importing input annotations: ")
    fusInfoDF = pd.read_csv(args.inFusInfo, sep="\t", skiprows=1,
                            names=['fusion_id','exon_id','gene_id','flag_multigene'],
                            usecols=['fusion_id','exon_id','gene_id'])
    fusBedDF = pd.read_csv(args.inFusBED, sep="\t", header=None,
                           names=['chr','fusion_start','fusion_stop','fusion_id','score','strand'],
                            usecols=['chr','fusion_start','fusion_stop','fusion_id','strand'])
    exonInfoDF = pd.read_csv(args.inExonInfo, sep=",", skiprows=1,
                             names=['chr','start','stop','strand','exon_id','transcript_id','gene_id'],
                             usecols=['chr','exon_id','transcript_id','gene_id'])
    # Add dataframes to SQL database
    cur.execute("create table fusInfoTable(fusion_id text, exon_id text, gene_id text);")
    cur.execute("create table fusBEDTable(chr text, fusion_start int, fusion_stop int, fusion_id text, strand text);")
    cur.execute("create table exonInfoTable (chr text, exon_id text, transcript_id text, gene_id text);")
    fusInfoDF.to_sql("fusInfoTable", con, if_exists="replace")
    fusBedDF.to_sql("fusBEDTable", con, if_exists="replace")
    exonInfoDF.to_sql("exonInfoTable", con, if_exists="replace")

    # Sort fusion tables and merge
    cur.execute("SELECT * FROM fusInfoTable ORDER BY fusion_id;")
    cur.execute("SELECT * FROM fusBEDTable ORDER BY fusion_id;")
    cur.execute("CREATE TABLE fusInfoTable2 AS SELECT fi.fusion_id, exon_id, gene_id, chr, "
                "fusion_start, fusion_stop, strand "
                "FROM fusInfoTable fi INNER JOIN fusBEDTable fb "
                "ON fi.fusion_id = fb.fusion_id ;")

    cur.execute("""Select * from exonInfoTable""")
    allEx2xs = cur.fetchall()
    cur.execute("""CREATE TABLE IF NOT EXISTS exon2xscript
                        (chr TEXT, exon_id TEXT, transcript_id TEXT, gene_id TEXT)""")
    for exon in allEx2xs:
        xscripts = exon[3].split("|")
        for xs in xscripts:
            cur.execute("INSERT INTO exon2xscript VALUES(:chr, :exon_id, :transcript_id, :gene_id)",
                           {"chr": exon[1], "exon_id": exon[2], "transcript_id": xs, "gene_id": exon[4]})

    # Merge exon2xsDF and fusInfoDF2
    cur.execute("SELECT * FROM fusInfoTable2 ORDER BY chr, gene_id, exon_id;")
    cur.execute("SELECT * FROM exon2xscript ORDER BY  chr, gene_id, exon_id;")
    cur.execute("CREATE TABLE fus2xscript AS SELECT fi.exon_id, fi.gene_id, fi.chr, fusion_id, "
                "fusion_start, fusion_stop, strand, transcript_id "
                "FROM fusInfoTable2 fi INNER JOIN exon2xscript ex "
                "ON fi.chr = ex.chr AND fi.gene_id = ex.gene_id AND fi.exon_id = ex.exon_id ;")

    print("Counting transcripts per gene to use for flagging unique, common, constitutive: " + str(datetime.datetime.now()))
    cur.execute("CREATE TABLE xs2gene AS SELECT distinct gene_id, transcript_id FROM exon2xscript ;")
    cur.execute("CREATE TABLE xsPerGene AS SELECT gene_id, count(gene_id) AS total_transcripts_per_genes FROM xs2gene  GROUP BY gene_id;")
    # Concatenate distinct exons, transcripts, and genes per fusion
    cur.execute("CREATE TABLE fusAnnot AS SELECT distinct fusion_id, chr, fusion_start, fusion_stop, strand, "
                "replace(group_concat(distinct exon_id),',', '|') AS exon_id, count(distinct exon_id) AS exons_per_fusion, "
                "replace(group_concat(distinct transcript_id),',', '|') AS transcript_id, count(distinct transcript_id) AS transcripts_per_fusion, "
                "replace(group_concat(distinct gene_id),',', '|') AS gene_id, count(distinct gene_id) AS genes_per_fusion "
                "FROM fus2xscript "
                "GROUP BY fusion_id ;")
    # Get max number of possible transcripts per gene(s) for each fusion and merge back into annotation
    cur.execute("CREATE TABLE fus2gene AS SELECT distinct fusion_id, gene_id FROM fus2xscript ;")
    cur.execute("SELECT * FROM fus2gene ORDER BY gene_id;")
    cur.execute("SELECT * FROM xsPerGene ORDER BY gene_id;")
    cur.execute("CREATE TABLE MaxXsPerFus AS SELECT fg.fusion_id AS fusion_id, sum(xg.total_transcripts_per_genes) AS total_transcripts_per_genes "
                "FROM fus2gene fg LEFT JOIN xsPerGene xg ON fg.gene_id = xg.gene_id GROUP BY fusion_id ;")
    cur.execute("SELECT * FROM MaxXsPerFus ORDER BY fusion_id;")
    cur.execute("SELECT * FROM fusAnnot ORDER BY fusion_id;")
    cur.execute("CREATE TABLE fusAnnot2 AS SELECT * "
                "FROM fusAnnot fa INNER JOIN MaxXsPerFus xf ON fa.fusion_id = xf.fusion_id ;")
    # Call annotation frequency, and flag if multigene
    cur.execute("ALTER TABLE fusAnnot2 ADD annotation_frequency text ;")
    cur.execute("ALTER TABLE fusAnnot2 ADD flag_multigene int ;")
    cur.execute("UPDATE fusAnnot2 "
                "SET annotation_frequency = (case when transcripts_per_fusion=0 then 'Unannotated' "
                "when transcripts_per_fusion=1 then 'Unique' "
                "when transcripts_per_fusion=total_transcripts_per_genes then 'Constitutive' "
                "when transcripts_per_fusion>1 and transcripts_per_fusion<total_transcripts_per_genes then 'Common' "
                "end);")
    cur.execute("UPDATE fusAnnot2 "
                "SET flag_multigene = (case when genes_per_fusion > 1 then '1'"
                "when genes_per_fusion = 1 then '0' end);")
    fusionAllInfoDF = pd.read_sql("SELECT * FROM fusAnnot2;", con)
    fusionAllInfoDF=fusionAllInfoDF[['fusion_id','chr', 'fusion_start','fusion_stop', 'strand','exons_per_fusion','exon_id',
                      'transcripts_per_fusion','transcript_id','genes_per_fusion','gene_id','annotation_frequency',
                      'flag_multigene','total_transcripts_per_genes']]
    fusionAllInfoBedDF=fusionAllInfoDF[['chr','fusion_start','fusion_stop','fusion_id']]
    # Write output
    with open(args.outCSV, 'w') as outCSVFile:
        fusionAllInfoDF.to_csv(outCSVFile, encoding='utf-8', index=False, sep=',')
    with open(args.outBED, 'w') as outBEDFile:
        fusionAllInfoBedDF.to_csv(outBEDFile, encoding='utf-8', index=False, sep='\t', header=None)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    main()
    logger.info('Script complete')
