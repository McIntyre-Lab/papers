#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2017-12-19
# NAME: import_and_format_exons.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script imports the junctions CSV file created by extractJunctions.py and reformats it for adding
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
    parser = argparse.ArgumentParser(description="Import exons CSV file")
    parser.add_argument("--input", dest="csvInput", action='store', required=True, help="Input CSV file")
    parser.add_argument("--output", dest="outCSV", action='store', required=True, help="Output CSV file name")

    args = parser.parse_args()
    return args

def main():
    # Import exon annotations and format
    logger.info("Importing CSV file of all annotated exons")
    exonDF = pd.read_csv(args.csvInput, sep=",")

    # Sort exons by chrom, gene, start and stop positions
    exonSortedDF = exonDF.sort_values(by=['chrom','gene_id','start','stop','exon_id'])

    # Create an index of chrom and gene_id to iterate over
    # Then use this to subset the exon dataframe
    geneIndex = exonSortedDF[['chrom', 'gene_id']]
    geneIndex = geneIndex.drop_duplicates(['chrom','gene_id'])
    ## Create sets of lists: one for exon-to-exon group, one for ref coordinates per group
    ## (i.e. start and stop of longest exon)
    # Leave this as pandas, since this would be too difficult to redo for SQLite3 (lots of table joining)
    # List set 1: gene_id, exon_group, refStart, refStop
    exonGeneList1=[]
    exonGroupList1=[]
    # List set 2: exon_id, gene_id, exon_group
    exonList=[]
    exonGeneList2=[]
    exonGroupList2=[]
    for indexGene, rowGene in geneIndex.iterrows():
        currentChrom = rowGene['chrom']
        currentGene = rowGene['gene_id']
        subsetExonDF=exonSortedDF.loc[(exonSortedDF['chrom'] == currentChrom) & (exonSortedDF['gene_id'] == currentGene)]
        exonGroupNum = 1
        exonIndexNum = 1
        for indexExon, rowExon in subsetExonDF.iterrows():
            # If first exon in gene
            if exonIndexNum == 1:
                exonGroupNum=exonGroupNum
                regionStart = rowExon['start']
                regionStop = rowExon['stop']
                exonLen= rowExon['stop'] - rowExon['start']
                currExonLen=exonLen
                exonID=rowExon['exon_id']
                exonIndexNum=exonIndexNum + 1
                # Append exon info to list
                exonList.append(exonID)
                exonGeneList2.append(currentGene)
                exonGroupList2.append(exonGroupNum)
            else : # If not the first exon in gene
                # Is next exon within the current region or in a new region?
                if rowExon['start'] <= regionStop : # Exon starts before region ends
                    if rowExon['stop'] >= regionStop : # Exon stops after region currently ends
                        regionStop=rowExon['stop']
                        regionStart=regionStart
                    else : # Exon stops before region currently ends
                        regionStop=regionStop
                        regionStart=regionStart
                    exonGroupNum = exonGroupNum
                    exonLen = rowExon['stop'] - rowExon['start']
                    # Append exon info to list
                    exonID = rowExon['exon_id']
                    exonList.append(exonID)
                    exonGeneList2.append(currentGene)
                    exonGroupList2.append(exonGroupNum)
                else : # Exon is in a new region
                    # First, write current exonGroup1 and exonGene1 info
                    exonGeneList1.append(currentGene)
                    exonGroupList1.append(exonGroupNum)
                    # Now get next set of information
                    exonGroupNum=exonGroupNum + 1
                    regionStart=rowExon['start']
                    regionStop=rowExon['stop']
                    exonID = rowExon['exon_id']
                    exonList.append(exonID)
                    exonGeneList2.append(currentGene)
                    exonGroupList2.append(exonGroupNum)
        # We have iterated through all exons in gene, so now write data to lists
        exonGeneList1.append(currentGene)
        exonGroupList1.append(exonGroupNum)
    # Open SQL database
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    # Put dataframes into database as tabless
    exon2GroupDF = pd.DataFrame({'exon_id': exonList, 'gene_id': exonGeneList2, 'exon_group': exonGroupList2})
    exonGroupsDF = pd.DataFrame({'gene_id': exonGeneList1, 'exon_group': exonGroupList1})
    cur.execute("CREATE TABLE IF NOT EXISTS exon2Group (exon_id TEXT, gene_id TEXT, exon_group TEXT); ")
    exon2GroupDF.to_sql("exon2Group", con, if_exists="replace")
    cur.execute("CREATE TABLE IF NOT EXISTS exonGroups (gene_id TEXT, exon_group TEXT); ")
    exonGroupsDF.to_sql("exonGroups", con, if_exists="replace")
    cur.execute("CREATE TABLE IF NOT EXISTS exonSorted (chrom TEXT, start INT, stop INT, strand TEXT,"
                "exon_id TEXT, transcript_id TEXT, gene_id TEXT); ")
    exonSortedDF.to_sql("exonSorted", con, if_exists="replace")
    # Sort and merge
    cur.execute("SELECT * FROM exon2Group ORDER BY gene_id, exon_group ; ")
    cur.execute("SELECT * FROM exonGroups ORDER BY gene_id, exon_group ; ")
    cur.execute("CREATE TABLE exon2Group2 AS SELECT in1.gene_id, in1.exon_group, in1.exon_id "
                "FROM exon2Group in1 INNER JOIN exonGroups in2 "
                "ON in1.gene_id = in2.gene_id AND in1.exon_group = in2.exon_group ;")
    cur.execute("SELECT * FROM exon2Group2 ORDER BY gene_id, exon_id ; ")
    cur.execute("SELECT * FROM exonSorted ORDER BY gene_id, exon_id ;")
    cur.execute("CREATE TABLE exon2Group3 AS SELECT in1.gene_id, in1.exon_id, in1.exon_group, in2.chrom, "
                "in2.start, in2.stop, in2.strand, in2.transcript_id "
                "FROM exon2Group2 in1 INNER JOIN exonSorted in2 "
                "ON in1.gene_id = in2.gene_id AND in1.exon_id = in2.exon_id ;")
    # Now identify exon groups that have alternative donors/acceptors, and which do not
    # Basically, count the number of distinct donor/acceptor sites per exon group
    cur.execute("CREATE TABLE site2group AS SELECT distinct gene_id, exon_group, start, stop FROM exon2Group3 ;")
    cur.execute("CREATE TABLE sitesPerGroup AS SELECT gene_id, exon_group, "
                "count(distinct stop) AS donors_per_group, count(distinct start) AS acceptors_per_group "
                "FROM site2group  GROUP BY gene_id, exon_group ; ")
    # Flag exon groups that have alternative donors or acceptors
    cur.execute("ALTER TABLE sitesPerGroup ADD flag_alt_donor INT ;")
    cur.execute("ALTER TABLE sitesPerGroup ADD flag_alt_acceptor INT ;")
    cur.execute("UPDATE sitesPerGroup "
                "SET flag_alt_donor = (case when donors_per_group > 1 then '1' "
                "when donors_per_group = 1 then '0' end);")
    cur.execute("UPDATE sitesPerGroup "
                "SET flag_alt_acceptor = (case when acceptors_per_group > 1 then '1' "
                "when acceptors_per_group = 1 then '0' end);")
    # Join with exon info, convert into pandas DF and export
    cur.execute("SELECT * FROM sitesPerGroup ORDER BY gene_id, exon_group ;")
    cur.execute("SELECT * FROM exon2Group3 ORDER BY gene_id, exon_group ;")
    cur.execute("CREATE TABLE exon2Group4 AS SELECT in1.*, in2.flag_alt_donor, in2.flag_alt_acceptor "
                "FROM exon2Group3 in1 INNER JOIN sitesPerGroup in2 "
                "ON in1.gene_id = in2.gene_id AND in1.exon_group = in2.exon_group ;")
    exon2GroupDF4 = pd.read_sql("SELECT * FROM exon2Group4;", con)
    # Export CSV of exons with donor/acceptor flags
    logger.info("Exporting formatted exon info")
    with open(args.outCSV, 'w') as outFile:
        exon2GroupDF4.to_csv(outFile, encoding='utf-8', index=False)

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
