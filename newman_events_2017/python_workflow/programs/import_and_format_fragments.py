#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-04
# NAME: import_and_format_fragments.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports the annotations for exon fragments and fusions (exonic regions) extracted from GFF
# and formats them in preparation for creating event-to-transcript index
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
    parser.add_argument("--input-fragment-file", dest="inFragInfo", action='store', required=True,
                        help="Input exon fragment CSV file")
    parser.add_argument("--input-fusion-file", dest="inFusInfo", action='store', required=True,
                        help="Input fusion info TSV file")
    parser.add_argument("--input-fusion-bed", dest="inFusBED", action='store', required=True,
                        help="Input fusion BED file")
    parser.add_argument("--input-exon-file", dest="inExonInfo", action='store', required=True,
                        help="Input exon info CSV file")
    parser.add_argument("--outCSV", dest="outCSV", action='store', required=True, help="Output CSV file name")
    parser.add_argument("--outBED", dest="outBED", action='store', required=True, help="Output BED file name")

    args = parser.parse_args()
    return args

def main():
    # Import exon annotations and format
    logger.info("Importing input annotations")
    fragInfoDF = pd.read_csv(args.inFragInfo, sep=",", skiprows=1,
                            names=['chr','region_start','region_stop','fragment_start','fragment_stop','exon_id'])
    fusInfoDF = pd.read_csv(args.inFusInfo, sep="\t", skiprows=1,
                            names=['fusion_id','exon_id','gene_id','flag_multigene'],
                            usecols=['fusion_id','exon_id'])
    fusBedDF = pd.read_csv(args.inFusBED, sep="\t", header=None,
                           names=['chr','fusion_start','fusion_stop','fusion_id','score','strand'],
                            usecols=['chr','fusion_start','fusion_stop','fusion_id'])
    exonInfoDF = pd.read_csv(args.inExonInfo, sep=",", skiprows=1,
                             names=['chr','start','stop','strand','exon_id','transcript_id','gene_id'],
                             usecols=['chr','exon_id','transcript_id','gene_id'])
    # Connect to database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Add dataframes to SQL database
    cur.execute("CREATE TABLE fragInfo (chr TEXT, region_start INT, region_stop INT, fragment_start INT, "
                "fragment_stop INT, exon_id TEXT) ;")
    fragInfoDF.to_sql("fragInfo", con, if_exists="replace")
    cur.execute("CREATE TABLE fusInfo (fusion_id TEXT, exon_id TEXT) ;")
    fusInfoDF.to_sql("fusInfo", con, if_exists="replace")
    cur.execute("CREATE TABLE fusBed (chr TEXT, fusion_start INT, fusion_stop INT, fusion_id TEXT ) ;")
    fusBedDF.to_sql("fusBed", con, if_exists="replace")
    cur.execute("CREATE TABLE exonInfo (chr TEXT, exon_id TEXT, transcript_id TEXT, gene_id TEXT) ;")
    exonInfoDF.to_sql("exonInfo", con, if_exists="replace")
    # Sort and merge fusion BED and info tables
    cur.execute("SELECT * FROM fusInfo ORDER BY fusion_id ; ")
    cur.execute("SELECT * FROM fusBed ORDER BY fusion_id ; ")
    cur.execute("CREATE TABLE fusInfo2 AS SELECT in1.*, in2.exon_id "
                "FROM fusBed in1 INNER JOIN fusInfo in2 "
                "ON in1.fusion_id = in2.fusion_id ; ")
    # Uncat exon_ids in fragments table
    cur.execute("Select * from fragInfo")
    allEx2Frag = cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS fragInfo2 (chr TEXT, exon_id TEXT, region_start INT, region_stop INT,"
                " fragment_start INT, fragment_stop INT ); ")
    for frag in allEx2Frag:
        exons = frag[6].split("|")
        for ex in exons:
            cur.execute("INSERT INTO fragInfo2 VALUES(:chr, :exon_id, :region_start, :region_stop, :fragment_start, "
                        ":fragment_stop )", {"chr": frag[1], "exon_id": ex, "region_start": frag[2],
                                             "region_stop": frag[3],"fragment_start": frag[4],
                                             "fragment_stop": frag[5] })
    # Join with fusion info on chr and exon
    cur.execute("SELECT * FROM fragInfo2 ORDER BY chr, exon_id, region_start, region_stop; ")
    cur.execute("SELECT * FROM fusInfo2 ORDER BY chr, exon_id, fusion_start, fusion_stop; ")
    cur.execute("CREATE TABLE frag2fusInfo AS SELECT in1.*, in2.fusion_id, in2.fusion_start, in2.fusion_stop "
                "FROM fragInfo2 in1 INNER JOIN fusInfo2 in2 "
                "ON in1.chr = in2.chr AND in1.exon_id = in2.exon_id AND in1.region_start = in2.fusion_start "
                "AND in1.region_stop = in2.fusion_stop ;")
    # Create a fragment ID based on the fusion ID
    cur.execute("CREATE TABLE fragCoord AS SELECT distinct fusion_id, fragment_start, fragment_stop "
                "FROM frag2fusInfo ORDER BY fusion_id, fragment_start, fragment_stop ; ")
    cur.execute("CREATE TABLE fragCoord2 AS SELECT in1.fusion_id as fusion_id, in1.fragment_start as fragment_start, "
                "in1.fragment_stop as fragment_stop, count(in2.fusion_id)+1 AS fragment_number "
                "FROM fragCoord in1 LEFT JOIN fragCoord in2 "
                "ON in1.fusion_id = in2.fusion_id AND in2.fragment_start < in1.fragment_start "
                "GROUP BY in1.fusion_id, in1.fragment_start ; ")
    cur.execute("CREATE TABLE fragCoord3 AS SELECT fusion_id, fragment_start, fragment_stop,"
                "fusion_id || ':' || fragment_number AS fragment_id "
                "FROM fragCoord2 ;")
    # Merge fragment IDs in with fragment info
    cur.execute("SELECT * FROM fragCoord3 ORDER BY fusion_id, fragment_start, fragment_stop ;")
    cur.execute("SELECT * FROM frag2fusInfo ORDER BY fusion_id, fragment_start, fragment_stop ;")
    cur.execute("CREATE TABLE frag2fusInfo2 AS SELECT in1.*, in2.fragment_id "
                "FROM frag2fusInfo in1 INNER JOIN fragCoord3 in2 "
                "ON in1.fusion_id = in2.fusion_id AND in1.fragment_start = in2.fragment_start "
                "AND in1.fragment_stop = in2.fragment_stop ; ")
    # Join exon info to get transcripts and genes
    cur.execute("SELECT * FROM frag2fusInfo2 ORDER BY chr, exon_id; ")
    cur.execute("SELECT * FROM exonInfo ORDER BY chr, exon_id; ")
    cur.execute("CREATE TABLE frag2exon AS SELECT in1.*, in2.transcript_id, in2.gene_id "
                "FROM frag2fusInfo2 in1 INNER JOIN exonInfo in2 "
                "ON in1.chr = in2.chr AND in1.exon_id = in2.exon_id ;")
    # Count the number of transcripts per gene. I will be using the to flag fusions as unique, common, constitutive
    # For each exon, uncat transcripts
    cur.execute("Select * from exonInfo")
    allxs2gene = cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS xs2gene (transcript_id TEXT, gene_id TEXT); ")
    for exon in allxs2gene:
        xscripts = exon[3].split("|")
        for xs in xscripts:
            cur.execute("INSERT INTO xs2gene VALUES(:transcript_id, :gene_id )",
                        {"transcript_id": xs, "gene_id": exon[4] })
    cur.execute("CREATE TABLE xsPerGene AS SELECT distinct gene_id, transcript_id FROM xs2gene ORDER BY gene_id, transcript_id ;")
    cur.execute("CREATE TABLE xsPerGene2 AS SELECT gene_id, count(gene_id) AS total_transcripts_per_genes FROM xsPerGene "
                "GROUP BY gene_id; ")
    # Count the number of distinct exons, transcripts and genes per fragment
    # First de-concat transcripts and genes
    cur.execute("SELECT * FROM frag2exon")
    allXs2Frag = cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS xs2frag (chr TEXT, exon_id TEXT, region_start INT, region_stop INT, fragment_start INT, "
                "fragment_stop INT, fusion_id TEXT, fusion_start INT, fusion_stop INT, fragment_id TEXT, transcript_id TEXT, "
                "gene_id TEXT); ")
    for frag in allXs2Frag:
        xscripts = frag[10].split("|")
        for xs in xscripts:
            cur.execute("INSERT INTO xs2frag VALUES(:chr, :exon_id, :region_start, :region_stop, :fragment_start, :fragment_stop, "
                        ":fusion_id, :fusion_start, :fusion_stop, :fragment_id, :transcript_id, :gene_id )",
                        {"chr": frag[0], "exon_id": frag[1], "region_start": frag[2], "region_stop": frag[3],
                         "fragment_start": frag[4], "fragment_stop": frag[5], "fusion_id": frag[6], "fusion_start": frag[7],
                         "fusion_stop": frag[8], "fragment_id": frag[9], "transcript_id": xs, "gene_id": frag[11]})
    cur.execute("SELECT * FROM xs2frag")
    allGene2Frag = cur.fetchall()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS gene2frag (chr TEXT, exon_id TEXT, region_start INT, region_stop INT, fragment_start INT, "
        "fragment_stop INT, fusion_id TEXT, fusion_start INT, fusion_stop INT, fragment_id TEXT, transcript_id TEXT, "
        "gene_id TEXT); ")
    for frag in allGene2Frag:
        genes = frag[11].split("|")
        for gn in genes:
            cur.execute(
                "INSERT INTO gene2frag VALUES(:chr, :exon_id, :region_start, :region_stop, :fragment_start, :fragment_stop, "
                ":fusion_id, :fusion_start, :fusion_stop, :fragment_id, :transcript_id, :gene_id )",
                {"chr": frag[0], "exon_id": frag[1], "region_start": frag[2], "region_stop": frag[3],
                 "fragment_start": frag[4], "fragment_stop": frag[5], "fusion_id": frag[6], "fusion_start": frag[7],
                 "fusion_stop": frag[8], "fragment_id": frag[9], "transcript_id": frag[10], "gene_id": gn})
    # Concatenate distinct exons, transcripts, and genes per fusion
    cur.execute("CREATE TABLE fragAnnot AS SELECT distinct fragment_id, chr, fragment_start, fragment_stop, fusion_id, "
                "fusion_start, fusion_stop, replace(group_concat(distinct exon_id),',', '|') AS exon_id, "
                "count(distinct exon_id) AS num_exons_per_frag, "
                "replace(group_concat(distinct transcript_id),',', '|') AS transcript_id, "
                "count(distinct transcript_id) AS num_xscripts_per_frag, "
                "replace(group_concat(distinct gene_id),',', '|') AS gene_id, "
                "count(distinct gene_id) AS num_genes_per_frag "
                "FROM gene2frag "
                "GROUP BY fragment_id ;")
    # Get max number of possible transcripts per gene(s) for each fusion and merge back into annotation
    cur.execute("CREATE TABLE frag2gene AS SELECT distinct fragment_id, gene_id FROM gene2frag ;")
    cur.execute("SELECT * FROM frag2gene ORDER BY gene_id;")
    cur.execute("SELECT * FROM xsPerGene2 ORDER BY gene_id;")
    cur.execute("CREATE TABLE MaxXsPerFrag AS SELECT fg.fragment_id AS fragment_id, sum(xg.total_transcripts_per_genes) AS total_transcripts_per_genes "
                "FROM frag2gene fg LEFT JOIN xsPerGene2 xg ON fg.gene_id = xg.gene_id GROUP BY fragment_id ;")
    cur.execute("SELECT * FROM MaxXsPerFrag ORDER BY fragment_id;")
    cur.execute("SELECT * FROM fragAnnot ORDER BY fragment_id;")
    cur.execute("CREATE TABLE fragAnnot2 AS SELECT * "
                "FROM fragAnnot fa INNER JOIN MaxXsPerFrag xf ON fa.fragment_id = xf.fragment_id ;")
    # Call annotation frequency, and flag if multigene
    cur.execute("ALTER TABLE fragAnnot2 ADD annotation_frequency text ;")
    cur.execute("ALTER TABLE fragAnnot2 ADD flag_multigene int ;")
    cur.execute("UPDATE fragAnnot2 "
                "SET annotation_frequency = (case when num_xscripts_per_frag=0 then 'Unannotated' "
                "when num_xscripts_per_frag=1 then 'Unique' "
                "when num_xscripts_per_frag=total_transcripts_per_genes then 'Constitutive' "
                "when num_xscripts_per_frag>1 and num_xscripts_per_frag<total_transcripts_per_genes then 'Common' "
                "end);")
    cur.execute("UPDATE fragAnnot2 "
                "SET flag_multigene = (case when num_genes_per_frag > 1 then '1'"
                "when num_genes_per_frag = 1 then '0' end);")
    fragAllInfoDF = pd.read_sql("SELECT * FROM fragAnnot2;", con)
    fragAllInfoDF = fragAllInfoDF[['fragment_id','chr','fragment_start','fragment_stop','fusion_id','fusion_start',
                                   'fusion_stop','exon_id','num_exons_per_frag','transcript_id','num_xscripts_per_frag',
                                   'gene_id','num_genes_per_frag','annotation_frequency','flag_multigene',
                                   'total_transcripts_per_genes']]
    # Export CSV of exon fragment info and a BED file for coverage
    logger.info("Exporting formatted exon fragment info")
    with open(args.outCSV, 'w') as outFile:
        fragAllInfoDF.to_csv(outFile, encoding='utf-8', index=False)
    logger.info("Exporting formatted exon fragment BED")
    fragAllInfoBedDF = fragAllInfoDF[['chr', 'fragment_start', 'fragment_stop', 'fragment_id']]
    with open(args.outBED, 'w') as outBEDFile:
        fragAllInfoBedDF.to_csv(outBEDFile, encoding='utf-8', index=False, sep='\t', header=None)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    main()
    logger.info('Script complete')
