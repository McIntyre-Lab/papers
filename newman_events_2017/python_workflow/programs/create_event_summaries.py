#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-04-16
# NAME: create_event_summaries.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script creates summaries for each event in the given input file. It takes a wide-formatted dataset of counts
# by event, annotations, detection flags, and a design file, and outputs a summary file, detailing group means,
# group detection, annotation frequency of events, transcripts and genes.
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
    parser.add_argument("-i", "--input-data", dest="inCounts", required=True, help="Wide-formatted dataset of counts")
    parser.add_argument("-d", "--design-file", dest="inDesign", required=True, help="Design file to relate samples to groups")
    parser.add_argument("-a", "--annotation-file", dest="inAnnot", required=True,
                        help="Formatted annotation file for events (make sure the correct annotation file is specified! ")
    parser.add_argument("-f", "--detection-flags-file", dest="inFlags", required=True, help="Detection flags for events")
    parser.add_argument("-j", "--junction-sequence-index", dest="inJuncSeq", required=False,
                        help="Junction-to-sequence index file created by extract_junction_sequence.py. This is REQUIRED "
                        "for junctions only. Exonic regions, fragments and introns do not need this infomation")
    # User-defined values
    parser.add_argument("-g", "--group-variable", dest="userGroup", required=True, 
                        help="Variable in design file used to group samples by treatment, condition, etc.")
    parser.add_argument("-l", "--event-length", dest="userLength", required=False, type=int, 
                        help="Minimum length (in bp) of event to be considered for analysis. Any event (exon fragment,"
                             "jnction, etc.) that is less than this value will be excluded from analysis. If not specified then "
                             "events are not filtered by their length.")

    # Outputs
    parser.add_argument("-o", "--output-summary", dest="outFlags", required=True, help="Output TSV of event detection flags by "
                                                                            "treatment group")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    # Import counts, design file, annotation and detection flags
    countsWideDF = pd.read_csv(args.inCounts, sep="\t")
    designDF = pd.read_csv(args.inDesign, sep="\t")
    annotDF = pd.read_csv(args.inAnnot, sep=",")
    flagsDF = pd.read_csv(args.inFlags, sep="\t")

    # For junctions, check that the user also suppose
    if 'junction_id' in annotDF.columns:
        if args.inJuncSeq : # check that a junction-to-sequence index was provided
            try:
                juncSeqDF = pd.read_csv(args.inJuncSeq, sep=",", usecols=['junction_id','sequence_id'])
                juncSeqDF.to_sql("junc2seqIndex", con, if_exists="replace")
            except IOError:
                print("Junction-to-sequence index file not specified or does not exist! Please specify this with the --junction-sequence-index option")
            except ValueError:
                print("Junction-to-sequence index file does not appear to be valid. Check your index file.")
            except :
                print("An unexpected error occurred while importing the junction-to-sequence index")
                raise
        else:
            print("Junction-to-sequence index file not specified! Please specify this with the --junction-sequence-index option")
            raise NameError("No junction-to-index file provided")
    # Convert to a tall dataset
    countsTallDF = pd.melt(countsWideDF, id_vars=['event_id'], var_name='sampleID', value_name='APN')
    # Send tall counts dataset and design file to SQL and merge
    countsTallDF.to_sql("countsTall", con, if_exists="replace")
    designDF.to_sql("designInfo", con, if_exists="replace")
    annotDF.to_sql("annotInfo", con, if_exists="replace")
    flagsDF.to_sql("flagsInfo", con, if_exists="replace")
    cur.execute("CREATE TABLE countsKey AS SELECT in1.*, in2.event_id, in2.APN "
                "FROM designInfo in1 INNER JOIN countsTall in2 "
                "ON in1.sampleID = in2.sampleID ;")
    # Calculate group means
    cur.execute("CREATE TABLE countsMean AS SELECT event_id, "
                + args.userGroup + ", avg(APN) as mean_apn FROM countsKey GROUP BY event_id, "+ args.userGroup + ";")
    # Put means side-by-side
    groupList = designDF[args.userGroup].drop_duplicates(keep='first').tolist()
    counter = 1
    for group in range(0,len(groupList)) :
        groupName = groupList[group]
        cur.execute("CREATE TABLE tempMean AS SELECT event_id, mean_apn AS mean_apn_"+groupName+" FROM countsMean "
                    "WHERE " + args.userGroup + " = ? ORDER BY event_id;", (groupName, ))
        if counter == 1:
            cur.execute("CREATE TABLE eventsMean AS SELECT * FROM tempMean ;")
            cur.execute("DROP TABLE tempMean;")
            counter=counter + 1
        else:
            cur.execute("CREATE TABLE eventsMean2 AS SELECT in1.*, in2.mean_apn_"
                        + groupName + " FROM eventsMean in1 INNER JOIN tempMean in2 "
                                      "ON in1.event_id = in2.event_id "
                                      "ORDER BY in1.event_id ;")
            cur.execute("DROP TABLE tempMean;")
            cur.execute("DROP TABLE eventsMean;")
            cur.execute("ALTER TABLE eventsMean2 RENAME TO eventsMean ;")
    # Merge in detection flags
    cur.execute("CREATE TABLE eventsMeansFlag AS SELECT * "
                "FROM eventsMean in1 INNER JOIN flagsInfo in2 ON in1.event_id = in2.event_id")
    # Create a list of columns to extract from annotations. We want (if possible) annotation frequency,
    # flag_multigene, transcript list, gene list, event_type
    annotList=['event_id']
    if 'annotation_frequency' in annotDF.columns:
        annotList.append('annotation_frequency')
    if 'flag_multigene' in annotDF.columns:
        annotList.append('flag_multigene')
    if 'transcript_id' in annotDF.columns:
        annotList.append('transcript_id')
    if 'gene_id' in annotDF.columns:
        annotList.append('gene_id')
    if 'flag_junction_annotated' in annotDF.columns:
        annotList.append('flag_junction_annotated')
    if 'flag_border_junction' in annotDF.columns:
        annotList.append('flag_border_junction')
    if 'flag_alt_donor' in annotDF.columns:
        annotList.append('flag_alt_donor')
    if 'flag_alt_acceptor' in annotDF.columns:
        annotList.append('flag_alt_acceptor')
    if 'flag_exonskip' in annotDF.columns:
        annotList.append('flag_exonskip')
    annotListStr = ', '.join(annotList)
    # Merge in annotations
    if 'fragment_id' in annotDF.columns:
        cur.execute("CREATE TABLE annotInfo2 AS SELECT *, fragment_id AS event_id, (fragment_stop-fragment_start) AS event_length FROM annotInfo ;" )
    elif 'intron_id' in annotDF.columns:
        cur.execute("CREATE TABLE annotInfo2 AS SELECT *, intron_id AS event_id, (intron_stop-intron_start) AS event_length FROM annotInfo ;")
    elif 'junction_id' in annotDF.columns:
        cur.execute("CREATE TABLE annotInfo2 AS SELECT *, junction_id AS event_id, (donor_stop-donor_start + acceptor_stop-acceptor_start) AS event_length FROM annotInfo ;")
    elif 'fusion_id' in annotDF.columns:
        cur.execute("CREATE TABLE annotInfo2 AS SELECT *, fusion_id AS event_id, (fusion_stop-fusion_start) AS event_length FROM annotInfo ;" )
    else:
        cur.execute("CREATE TABLE annotInfo2 AS SELECT * FROM annotInfo ;")
    if args.userLength:
        cur.execute("CREATE TABLE annotInfo3 AS SELECT " + annotListStr + " FROM annotInfo2 WHERE event_length >= ? ORDER BY event_id ;", (args.userLength, ) )
    else :
        cur.execute("CREATE TABLE annotInfo3 AS SELECT " + annotListStr + " FROM annotInfo2 ORDER BY event_id ;")
    if 'junction_id' in annotDF.columns: # If junctions, we need to convert these from sequnce IDs to junction IDs
        cur.execute("CREATE TABLE eventsMeansFlag2 AS SELECT *, event_id AS sequence_id FROM eventsMeansFlag ;")
        cur.execute("CREATE TABLE eventsMeansFlag3 AS SELECT in1.*, in2.junction_id "
                    "FROM eventsMeansFlag2 in1 INNER JOIN junc2seqIndex in2 "
                    "ON in1.sequence_id = in2.sequence_id")
        cur.execute("CREATE TABLE eventsMeansFlagAnnot AS SELECT in1.*, in2.* "
                    "FROM annotInfo3 in1 INNER JOIN eventsMeansFlag3 in2 "
                    "ON in1.event_id = in2.junction_id ;")
        eventSummaryDF = pd.read_sql("SELECT * FROM eventsMeansFlagAnnot;", con)
        eventSummaryDF = eventSummaryDF.drop(['event_id:1','event_id:2','sequence_id','junction_id','index'], axis=1)
    else:
        cur.execute("CREATE TABLE eventsMeansFlagAnnot AS SELECT in1.*, in2.* "
                    "FROM eventsMeansFlag in1 INNER JOIN annotInfo3 in2 "
                    "ON in1.event_id = in2.event_id ;")
        eventSummaryDF = pd.read_sql("SELECT * FROM eventsMeansFlagAnnot;", con)
        eventSummaryDF = eventSummaryDF.drop(['event_id:1','event_id:2','index'], axis=1)
    # Write output flags file
    with open(args.outFlags, 'w') as outFile:
        eventSummaryDF.to_csv(outFile, encoding='utf-8', index=False, sep="\t")

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

