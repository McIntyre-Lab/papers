#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-04-16
# NAME: create_event_summaries.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script creates an output table that flags whether a gene is expressed or not, and whether it has
# any associated exonic regions of genic ambiguity (multigene fusions). The inputs are a design file and the fusion
# summary generated with create_event_summaries.py. The user must also specify which variable in there design file is
# to be used to group sample. The output is a summary table that includes a flag for each group in your design file
# indicating if a gene demonstrates evidence of expression, and a single flag that indicates whether the gene has at
# least one exonic region (fusion) of genic ambiguity This output table is used when creating the final transcript-
# level summary and intron-border summary, but excluding transcripts of these genes from the summary outputs
#
# REQUIRED PACKAGES: pandas    (tested with v0.19.2)
#                    argparse  (tested with v1.1)
#                    logging   (tested with v0.5.1.2)
#                    sqlite3
#
#######################################################################################################################

# Import required packages
import pandas as pd
import logging
import sqlite3
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes wide-formatted counts matrix and a design file and outputs a "
                                                 "TSV of group-wise detection flags")
    # Inputs
    parser.add_argument("-i", "--input-exonic-region-summary", dest="inSummary", required=True, help="Summary of fusions/exonic regions as genereated by create_event_summaries.py")
    parser.add_argument("-d", "--design-file", dest="inDesign", required=True, help="Design file to relate samples to groups")
    # User-defined values
    parser.add_argument("-g", "--group-variable", dest="userGroup", required=True, help="Variable in design file used to group samples by treatment, condition, etc.")

    # Outputs
    parser.add_argument("-o", "--output-gene-summary", dest="outSummary", required=True, help="Output TSV of gene-level expression flags")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()
    # Import fusion summary and design file
    fusionDF = pd.read_csv(args.inSummary, sep="\t")
    designDF = pd.read_csv(args.inDesign, sep="\t")
    # Send fusion summary DF to SQL database
    fusionDF.to_sql("fusionInfo", con, if_exists="replace")

    # First de-concatenate genes
    cur.execute("SELECT distinct event_id, gene_id, flag_multigene FROM fusionInfo ORDER BY event_id ; ")
    fus2gene= cur.fetchall()
    cur.execute("CREATE TABLE IF NOT EXISTS fus2gene (event_id TEXT, gene_id TEXT, flag_multigene INT);")
    for fusion in fus2gene:
        genes = fusion[1].split("|")
        for gn in genes:
            cur.execute("INSERT INTO fus2gene VALUES(:event_id, :gene_id, :flag_multigene ) ;",
                        {"event_id": fusion[0], "gene_id": gn, "flag_multigene": fusion[2]})
    # Get the list of genes that have at least one multigene region
    cur.execute("CREATE TABLE genesFlagMulti AS SELECT distinct gene_id, flag_multigene "
                "FROM fus2gene WHERE flag_multigene = 1 ORDER BY gene_id; ")
    # For each treatment group, count the number of detected fusions per gene
    # First merge fus2gene and
    # Create a list of groups to iterate over
    groupList = designDF[args.userGroup].drop_duplicates(keep='first').tolist()
    # Iterate over groups
    counter = 1
    for group in range(0,len(groupList)) :
        groupName = groupList[group]
        # Get detection flags for treatment group
        cur.execute("CREATE TABLE tempFusDtct AS SELECT event_id, flag_"
                    +groupName+"_detected FROM fusionInfo ORDER BY event_id ;")
        # Join with fus2gene
        cur.execute("CREATE TABLE tempFusDtct2gene AS SELECT in1.*, in2.gene_id "
                    "FROM tempFusDtct in1 INNER JOIN fus2gene in2 "
                    "ON in1.event_id = in2.event_id ;")
        # Count the number of detected fusions in each group
        cur.execute("CREATE TABLE tempFusCount AS SELECT distinct gene_id, count(distinct event_id) AS "
                    "num_exonic_regions, sum(flag_"+groupName+"_detected) AS num_exonic_regions_detected_"+groupName+" "
                    "FROM tempFusDtct2gene GROUP BY gene_id ;" )
        # Flag if a gene has at least one detected fusion
        cur.execute("ALTER TABLE tempFusCount ADD flag_gene_expressed_"+groupName+" INT ;")
        cur.execute("UPDATE tempFusCount SET flag_gene_expressed_"
                    +groupName+" = (CASE WHEN num_exonic_regions_detected_"+groupName+" > 0 THEN '1' "
                    "ELSE '0' END);")
        # Append
        if counter == 1:
            cur.execute("CREATE TABLE geneExp AS SELECT * FROM tempFusCount ;")
            cur.execute("DROP TABLE tempFusCount;")
            cur.execute("DROP TABLE tempFusDtct;")
            cur.execute("DROP TABLE tempFusDtct2gene;")
            counter=counter + 1
        else:
            cur.execute("CREATE TABLE geneExp2 AS SELECT in1.*, in2.num_exonic_regions_detected_"
                        + groupName + ", in2.flag_gene_expressed_"+groupName+" "
                        "FROM geneExp in1 INNER JOIN tempFusCount in2 "
                        "ON in1.gene_id = in2.gene_id "
                        "ORDER BY in1.gene_id ;")
            cur.execute("DROP TABLE geneExp;")
            cur.execute("DROP TABLE tempFusCount;")
            cur.execute("DROP TABLE tempFusDtct;")
            cur.execute("DROP TABLE tempFusDtct2gene;")
            cur.execute("ALTER TABLE geneExp2 RENAME TO geneExp ;")
    # Merge with the list of genes with at least one multigene fusion
    cur.execute("CREATE TABLE geneExpSummary AS SELECT in1.*, in2.flag_multigene AS flag_gene_has_multigene_exon "
                "FROM geneExp in1 LEFT JOIN genesFlagMulti in2 "
                "ON in1.gene_id = in2.gene_id; ")
    cur.execute("UPDATE geneExpSummary SET flag_gene_has_multigene_exon = 0 "
                "WHERE coalesce(flag_gene_has_multigene_exon, '') = '' ;")
    # Output summary file
    geneFlagsDF = pd.read_sql("SELECT * FROM geneExpSummary;", con)
    with open(args.outSummary, 'w') as outFile:
        geneFlagsDF.to_csv(outFile, encoding='utf-8', index=False, sep="\t")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    # Setting up logger
    logger = logging.getLogger()
    logger.info('Starting script')
    # Calling main script
    main()
    logger.info('Script complete')

