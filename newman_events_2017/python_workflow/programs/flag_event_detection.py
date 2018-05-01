#!/usr/bin/env python3
#######################################################################################################################
#
# DATE: 2018-04-13
# NAME: flag_event_detection.py
# AUTHOR: Jeremy R. B. Newman (jrbnewman@ufl.edu)
#
# DESCRIPTION: This script flags events (exons, exon fragments, introns, junctions, etc.) that are detected given some
# simple user-defined rules. It takes a wide-formatted matrix of counts (rows as events, columns as samples) and a
# design file relating columns to sample metadata such as treatment group. The user must specify what groups are
# present in their data, such that events can be flagged as detected by group. If there is no group information
# (i.e. user wants to consider all samples together), then a dummy-group variable must be provided.
# samples are treated as belonging to the same group. The user can also specify the minimum acceptable abundance value
# (e.g. APN, average depth per nucleotide) for an event to be considered "detected" in each group, as well as the minimum
# proportion of samples per group where that event is detected, for the event to be considered "detected" in the
# entire group. For example, if a user specifies a minimum abundance of 5, and a minimum proportion of 50%, this means
# that at least 50% of samples within a group must have an abundance of at least 5 for that event to be considered
# detected in that group. If "0" is supplied for either parameter, this will be interpreted as "> 0", otherwise this
# will be interpreted as ">= X", where X is the parameter value. If no parameters are specified then the default
# value of abundance cut-off is 0 (APN>0) and a default value of 0.5 for proportion of samples in which the event is
# detected.
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
    # Input data
    parser.add_argument("-i", "--input-data", dest="inCounts", required=True, help="Wide-formatted dataset of counts")
    parser.add_argument("-d", "--design-file", dest="inDesign", required=True, help="Design file to relate samples to groups")

    # User-defined parameters
    parser.add_argument("-g", "--group-variable", dest="userGroup", required=True, help="Variable name in design file to group samples on")
    parser.add_argument("-a", "--minimum-abundance", dest="userAPN", required=False,
                        help="Abundance must be greater than or equal to this value to call an event detected. If 0, then abundance must be "
                             "greater than this value")
    parser.add_argument("-p", "--minimum-proportion-detected", dest="userProp", required=False,
                        help="Minimum proportion (scale of 0 to 1) of samples in which the event is detected for the event to considered "
                             "detected in that group. Default is 0.5 (i.e. 50%)")

    # Output TSV
    parser.add_argument("-o", "--output-flags", dest="outFlags", required=True, help="Output TSV of event detection flags by "
                                                                            "treatment group")
    args = parser.parse_args()
    return args

def main():
    # Connect to SQL database
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    # Import counts data and design file
    countsWideDF = pd.read_csv(args.inCounts, sep="\t")
    designDF = pd.read_csv(args.inDesign, sep="\t")
    # Convert to a tall dataset
    countsTallDF = pd.melt(countsWideDF, id_vars=['event_id'], var_name='sampleID', value_name='APN')
    # Send tall counts dataset and design file to SQL and merge
    countsTallDF.to_sql("countsTall", con, if_exists="replace")
    designDF.to_sql("designInfo", con, if_exists="replace")
    cur.execute("CREATE TABLE countsKey AS SELECT in1.*, in2.event_id, in2.APN "
                "FROM designInfo in1 INNER JOIN countsTall in2 "
                "ON in1.sampleID = in2.sampleID ;")
    # Flag by event and sample, based on user-defined (or default) parameters
    ## Get user-defined thresholds. If they don't exist then set to 0
    if args.userAPN :
        minAPN = float(args.userAPN)
    else :
        minAPN = 0
    if args.userProp :
        minProp = float(args.userProp)
    else :
        minProp = 0
    # First, flag each individual observation per sample if above the user-defined APN threshold
    cur.execute("ALTER TABLE countsKey ADD flag_detected INT;")
    if minAPN == 0:
        cur.execute("UPDATE countsKey SET flag_detected = (CASE WHEN APN > 0 THEN '1' ELSE '0' END);")
    else:
        cur.execute("UPDATE countsKey SET flag_detected = (CASE WHEN APN >= ? THEN '1' ELSE '0' END);", (minAPN,))
    # Summarize detection per group
    cur.execute("CREATE TABLE countsFlagged AS SELECT event_id, "+args.userGroup+", avg(flag_detected) AS perc_trt FROM countsKey GROUP BY event_id, "+args.userGroup+" ;")
    # Create a list of groups to iterate over
    groupList = designDF[args.userGroup].drop_duplicates(keep='first').tolist()
    counter = 1
    for group in range(0,len(groupList)) :
        groupName = groupList[group]
        cur.execute("CREATE TABLE flagTemp AS SELECT event_id, perc_trt FROM countsFlagged "
                    "WHERE " + args.userGroup + " = ? ORDER BY event_id;", (groupName, ))
        cur.execute("ALTER TABLE flagTemp ADD flag_" + groupName + "_detected INT ;")
        if minProp == 0:  # min Prop = 0
            cur.execute("UPDATE flagTemp SET flag_" + groupName + "_detected = (CASE WHEN perc_trt > 0 THEN '1' ELSE '0' END) ;")
        else:
            cur.execute("UPDATE flagTemp SET flag_" + groupName +
                            "_detected = (CASE WHEN perc_trt >= ? THEN '1' ELSE '0' END) ;", (minProp, ))
        if counter == 1: # if first group value, then we create the output detection flags summary with group 1
            cur.execute("CREATE TABLE flagDtct AS SELECT event_id, flag_" + groupName + "_detected FROM flagTemp ; ")
            cur.execute("DROP TABLE flagTemp;")
            counter = counter + 1
        else:   # if not group 1, merge flagTemp with flags summary
            cur.execute("CREATE TABLE flagDtct2 AS SELECT in1.*, in2.flag_"
                        + groupName + "_detected FROM flagDtct in1 INNER JOIN flagTemp in2 "
                                      "ON in1.event_id = in2.event_id "
                                      "ORDER BY in1.event_id ;")
            cur.execute("DROP TABLE flagDtct;")
            cur.execute("DROP TABLE flagTemp;")
            cur.execute("ALTER TABLE flagDtct2 RENAME TO flagDtct ;")

    # Write output flags file
    detectFlagDF = pd.read_sql("SELECT * FROM flagDtct;", con)
    
    with open(args.outFlags, 'w') as outFile:
        detectFlagDF.to_csv(outFile, encoding='utf-8', index=False, sep="\t")

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

