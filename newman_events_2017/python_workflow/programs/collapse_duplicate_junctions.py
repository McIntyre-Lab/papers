#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-01-04
# NAME: collapse_duplicate_junctions.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script takes the complete set of logical junctions with donor and acceptor exon information and
# and other annotations and collapses it to its genomic coordinate, in the format chr:donor_site:acceptor_site:strand
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
from openpyxl.descriptors.base import ASCII


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import fragments and fusion annotation files")
    parser.add_argument("--input-junction-file", dest="inJuncInfo", action='store', required=True,
                        help="Input junction CSV file with donor/acceptor information")
    parser.add_argument("--input-exon-annotation", dest="inExonInfo", action='store', required=True,
                        help="Input exon annotation CSV, used for determining annotation frequency of junctions"
                             "(i.e. unique, common, or constitutive")
    parser.add_argument("--output-collapsed-junctions", dest="outCSV", action='store', required=True,
                        help="Output CSV file name")

    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Import uncollapsed junction annotations
    juncInDF = pd.read_csv(args.inJuncInfo, sep=",", skiprows=1,
                           names=['event_id','event_type','chr','strand','num_transcripts','transcript_id',
                                  'donor_exon','acceptor_exon','donor_start','donor_stop','donor_size',
                                  'acceptor_start','acceptor_stop','acceptor_size','flag_junction_annotated',
                                  'flag_border_junction','flag_exonskip','num_skipped_exons','cat_skipped_exons',
                                  'donor_gene','donor_group','flag_alt_donor','acceptor_gene','acceptor_group',
                                  'flag_alt_acceptor'], dtype={'flag_junction_annotated':int,'flag_border_junction':int,
                                                               'flag_exonskip':int, 'flag_alt_donor': int, 'flag_alt_acceptor':int})
    xsInfoDF = pd.read_csv(args.inExonInfo, sep=",", skiprows=1,
                             names=['chr','start','stop','strand','exon_id','transcript_id','gene_id'],
                             usecols=['transcript_id','gene_id'])
    # Convert dataframes to SQL tables
    juncInDF.to_sql("juncInfo", con, if_exists="replace")
    xsInfoDF.to_sql("xsInfo", con, if_exists="replace")

    # Count the number of distinct annotated transcripts per gene. Will need this for determining whether a junction is
    # unique, common, constitutive, or unannotated
    # Uncat and take distinct transcripts, genes
    cur.execute("Select distinct * from xsInfo")
    allXs2Gene = cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS xs2gene (transcript_id TEXT, gene_id TEXT);")
    for junc in allXs2Gene:
        xscripts = junc[1].split("|")
        for xs in xscripts:
            cur.execute("INSERT INTO xs2gene VALUES(:transcript_id, :gene_id) ;",
                           {"transcript_id": xs, "gene_id": junc[2]})
    cur.execute("CREATE TABLE xs2gene2 AS SELECT distinct transcript_id, gene_id FROM xs2gene ; ")
    cur.execute("CREATE TABLE xsPerGene AS SELECT gene_id, count(gene_id) AS total_transcripts_per_gene "
                "FROM xs2gene2 GROUP BY gene_id;")
    # Create a unique junction ID with the format: chr:donor:acceptor:strand, and add a temp "is_junction" flag
    cur.execute("CREATE TABLE juncInfo2 AS SELECT *, chr || ':' || donor_stop || ':' || acceptor_start || ':' || strand AS junction_id "
                "FROM juncInfo ORDER BY chr, donor_stop, acceptor_start, strand ;")
    cur.execute("ALTER TABLE juncInfo2 ADD flag_is_junction INT;")
    cur.execute("UPDATE juncInfo2 SET flag_is_junction = (CASE WHEN flag_border_junction=0 THEN '1' "
                "WHEN flag_border_junction=1 THEN '0' end); ")
    # For each distinct junction, take the minimum donor_start site, and maximum acceptor_stop site
    cur.execute("CREATE TABLE juncCoord AS SELECT junction_id, chr, min(donor_start) as donor_start, donor_stop, "
                "acceptor_start, max(acceptor_stop) as acceptor_stop, strand "
                "FROM juncInfo2 GROUP BY junction_id ;")
    # Concatenate and count the number of distinct transcripts, genes, donors, acceptors.
    # Set "unannotated" in transcript_id, and "intron" in donor/acceptor exon to empty for now (add this back in later)
    cur.execute("CREATE TABLE junc2ex2xs2gene AS SELECT junction_id, transcript_id, donor_exon, acceptor_exon, "
                "donor_gene, acceptor_gene FROM juncInfo2 ORDER BY junction_id ;")
    # Uncat transcripts
    cur.execute("Select * from junc2ex2xs2gene")
    allXs2Junc= cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS junc2ex2xs2gene2 (junction_id TEXT, transcript_id TEXT, donor_exon TEXT, "
                "acceptor_exon TEXT, donor_gene TEXT, acceptor_gene  TEXT);")
    for junc in allXs2Junc:
        xscripts = junc[1].split("|")
        for xs in xscripts:
            cur.execute("INSERT INTO junc2ex2xs2gene2 VALUES(:junction_id, :transcript_id, :donor_exon, :acceptor_exon, "
                        " :donor_gene, :acceptor_gene ) ;",
                        {"junction_id": junc[0], "transcript_id": xs, "donor_exon":junc[2], "acceptor_exon": junc[3],
                         "donor_gene":junc[4], "acceptor_gene":junc[5] })
    # Add gene_id from either donor or acceptor exon
    cur.execute("ALTER TABLE junc2ex2xs2gene2 ADD gene_id TEXT ;")
    cur.execute("UPDATE junc2ex2xs2gene2 "
                "SET gene_id = (CASE WHEN donor_exon='intron' THEN acceptor_gene "
                "WHEN donor_exon='intron' THEN acceptor_gene ELSE donor_gene END ) ; ")
    # Replace "unannotated" and "intron" from transcript_id/exon IDs to NULL
    cur.execute("UPDATE junc2ex2xs2gene2 SET donor_exon = NULL WHERE donor_exon='intron'; ")
    cur.execute("UPDATE junc2ex2xs2gene2 SET acceptor_exon = NULL WHERE acceptor_exon ='intron'; ")
    cur.execute("UPDATE junc2ex2xs2gene2 SET transcript_id = NULL WHERE transcript_id='Unannotated'; ")
    # Concatenate distinct donors, acceptors, transcripts and genes
    cur.execute(
        "CREATE TABLE juncInfoCollapsed AS SELECT distinct junction_id, replace(group_concat(distinct donor_exon),',', '|') as donor_exon_id, "
        "count(distinct donor_exon) AS num_donor_exons, replace(group_concat(distinct acceptor_exon),',', '|') as acceptor_exon_id, "
        "count(distinct acceptor_exon) AS num_acceptor_exons, replace(group_concat(distinct transcript_id),',', '|') as transcript_id, "
        "count(distinct transcript_id) AS num_transcripts, replace(group_concat(distinct gene_id),',', '|') as gene_id, "
        "count(distinct gene_id) AS num_genes FROM junc2ex2xs2gene2 GROUP BY junction_id ;")
    #(5) Count the total possible number of genes (i.e. sum total_xs per distinct gene)
    cur.execute("CREATE TABLE junc2gene AS SELECT distinct junction_id, gene_id FROM junc2ex2xs2gene2 ORDER BY gene_id; ")
    cur.execute("CREATE TABLE junc2XsPerGene AS SELECT in1.*, in2.total_transcripts_per_gene "
                "FROM junc2gene in1 INNER JOIN xsPerGene in2 ON in1.gene_id = in2.gene_id ;")
    cur.execute("CREATE TABLE juncMaxXs AS SELECT distinct junction_id, sum(total_transcripts_per_gene) AS "
                "total_transcripts_per_gene FROM junc2XsPerGene GROUP BY junction_id; ")
    #(6) Concatenate and count the number of distinct skipped exons
    cur.execute("Select distinct junction_id, cat_skipped_exons FROM juncInfo2 "
                "WHERE cat_skipped_exons != '(none)' ORDER BY junction_id ; ")
    allskip2Junc= cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS junc2skip (junction_id TEXT, exon_id TEXT);")
    for junc in allskip2Junc:
        exons = junc[1].split("|")
        for ex in exons:
            cur.execute("INSERT INTO junc2skip VALUES(:junction_id, :exon_id ) ;",
                        {"junction_id": junc[0], "exon_id": ex })
    cur.execute("CREATE TABLE juncCatSkip AS SELECT distinct junction_id, replace(group_concat(distinct exon_id),',', '|') "
        "AS skipped_exon_id, count(distinct exon_id) AS num_skipped_exons "
        "FROM junc2skip GROUP BY junction_id ;")
    # Collapse flags
    cur.execute("CREATE TABLE juncFlags AS SELECT junction_id, max(flag_junction_annotated) AS flag_junction_annotated, "
                "max(flag_is_junction) AS flag_is_junction, max(flag_border_junction) AS flag_border_junction, "
                "max(flag_exonskip) AS flag_exonskip, max(flag_alt_donor) AS flag_alt_donor, "
                "max(flag_alt_acceptor) AS flag_alt_acceptor FROM juncInfo2 GROUP BY junction_id;")
    # If junction is both flagged as an exon junction and as a border junction, then reclassify it as
    # an exon junction ONLY, as border junctions are supposed to have unambiguous intron sequence
    # Then also add back in the event type (exon_junction, exon_intron_border)
    cur.execute("UPDATE juncFlags SET flag_border_junction = 0 WHERE flag_is_junction = 1; ")
    cur.execute("ALTER TABLE juncFlags ADD event_type TEXT; ")
    cur.execute("UPDATE juncFlags SET event_type = (CASE WHEN flag_border_junction=1 THEN 'exon_intron_border' "
                "ELSE 'exon_junction' END) ;")
    # Merge it all together, then call junctions as unique, common, constitutive, unannotated
    cur.execute("CREATE TABLE juncAnnotWInfo AS SELECT in1.*, in2.donor_exon_id, in2.num_donor_exons, "
                "in2.acceptor_exon_id, in2.num_acceptor_exons, in2.transcript_id, in2.num_transcripts, "
                "in2.gene_id, in2.num_genes FROM juncCoord in1 INNER JOIN juncInfoCollapsed in2 "
                "ON in1.junction_id = in2.junction_id ORDER BY in1.junction_id ;")
    cur.execute("CREATE TABLE juncAnnotWMaxXs AS SELECT in1.*, in2.total_transcripts_per_gene "
                "FROM juncAnnotWInfo in1 INNER JOIN juncMaxXs in2 "
                "ON in1.junction_id = in2.junction_id ORDER BY in1.junction_id ;")
    cur.execute("CREATE TABLE juncAnnotWFlags AS SELECT in1.*, in2.flag_junction_annotated, in2.flag_border_junction, "
                "in2.flag_exonskip, in2.flag_alt_donor, in2.flag_alt_acceptor "
                "FROM juncAnnotWMaxXs in1 INNER JOIN juncFlags in2 "
                "ON in1.junction_id = in2.junction_id ORDER BY in1.junction_id ;")
    cur.execute("CREATE TABLE juncAnnotAll AS SELECT in1.*, in2.skipped_exon_id, in2.num_skipped_exons "
                "FROM juncAnnotWFlags in1 LEFT JOIN juncCatSkip in2 "
                "ON in1.junction_id = in2.junction_id ORDER BY in1.junction_id ; ")
    # Replace missing values: transcript_id "Unannotated", donor_exon_id "intron", acceptor_exon_id "intron"
    # skipped_exon_id "", num_skipped_exons "0", num_donor_exons "0", num_acceptor_exons "0"
    cur.execute("UPDATE juncAnnotAll SET transcript_id = 'Unannotated' WHERE coalesce(transcript_id, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET donor_exon_id ='intron' WHERE coalesce(donor_exon_id, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET acceptor_exon_id = 'intron' WHERE coalesce(acceptor_exon_id, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET skipped_exon_id = '' WHERE coalesce(skipped_exon_id, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET num_skipped_exons='0' WHERE coalesce(num_skipped_exons, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET num_donor_exons='0' WHERE coalesce(num_donor_exons, '') = ''; ")
    cur.execute("UPDATE juncAnnotAll SET num_acceptor_exons='0' WHERE coalesce(num_acceptor_exons, '') = ''; ")
    # Add in annotation frequency (U/C/X/N) and flag_multigene
    cur.execute("ALTER TABLE juncAnnotAll ADD annotation_frequency TEXT; ")
    cur.execute("ALTER TABLE juncAnnotAll ADD flag_multigene INT; ")
    cur.execute("UPDATE juncAnnotAll SET annotation_frequency = (CASE WHEN num_transcripts = 0 THEN 'Unannotated' "
                "WHEN num_transcripts = 1 THEN 'Unique' "
                "WHEN num_transcripts = total_transcripts_per_gene THEN 'Constitutive' "
                "ELSE 'Common' END); ")
    cur.execute("UPDATE juncAnnotAll SET flag_multigene = (CASE WHEN num_genes > 1 THEN '1' ELSE '0' END); ")
    # Export
    juncAllDF = pd.read_sql("SELECT junction_id, chr, donor_start, donor_stop, acceptor_start, acceptor_stop, "
                            "strand, num_transcripts, transcript_id, num_genes, gene_id, total_transcripts_per_gene, "
                            "annotation_frequency, flag_multigene, num_donor_exons, donor_exon_id, num_acceptor_exons, "
                            "acceptor_exon_id, flag_junction_annotated, flag_border_junction, flag_alt_donor, "
                            "flag_alt_acceptor, flag_exonskip, skipped_exon_id, num_skipped_exons "
                            "FROM juncAnnotAll ; ", con)
    juncAllDF = juncAllDF[['junction_id','chr','donor_start','donor_stop','acceptor_start','acceptor_stop','strand','num_transcripts',
         'transcript_id','num_genes','gene_id','total_transcripts_per_gene','annotation_frequency','flag_multigene',
         'num_donor_exons','donor_exon_id','num_acceptor_exons','acceptor_exon_id','flag_junction_annotated',
         'flag_border_junction','flag_alt_donor','flag_alt_acceptor','flag_exonskip','skipped_exon_id',
         'num_skipped_exons']]
    logger.info("Exporting collapsed junction annotation CSV")
    with open(args.outCSV, 'w') as outFile:
        juncAllDF.to_csv(outFile, encoding='utf-8', index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    logger = logging.getLogger()
    logger.info('Starting script')
    main()
    logger.info('Script complete')