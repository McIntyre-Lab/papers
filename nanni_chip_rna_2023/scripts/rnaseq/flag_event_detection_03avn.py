#!/usr/bin/env python3

### MODIFIED VERSION INFO
# The following arguments were added:
#     -f, --additional-flags which takes in a list of alternative groups (based on design file columns names) to make
#         flags for, which will have the same minimum value and proportion applied across all individual samples within
#         the groups (e.g. for 2 genotypes in a species, 2 sexes, 2 treatments, and 3 replicates = 12 samples per sex,
#         with a proportion of >50% with APN>0 then the detection of a species_sex group flag would require at least
#         7 of the 12 samples must have APN>0 to be detected)
#         *** Code for this is just a copy of the code used for the initial group (-g) put within a for loop and merged
#
#     -s, --all-sampleID-flags which takes in a list of groups (based on design file columns names) to make
#         flags for, which will require that all replicates to be on (aka proportion = 100%) or off (0%) with the 
#         minimum value (e.g. APN>0)
#         *** Code for this is just a copy of the code used for the initial group (-g) with proportion set to 1 or 0 
#             and put within a for loop and merged
#


# ORIGINAL SCRIPT INFO
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
    parser.add_argument("-f", "--additional-flags", dest="userFlags", required=False, action='append',
                        help="Additional flags based on the initial set of flags can be specified by giving the column name in the design "
                             "file for the new group to be flagged (will use the same abundance and proportion)")
    parser.add_argument("-s", "--all-sampleID-flags", dest="userAll", required=False, action='append',
                        help="Additional flags based on all sampleID's present in a group can be specified by giving the column name in the design "
                             "file for the new group to be flagged")

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
    countsWideDF = pd.read_csv(args.inCounts)
    designDF = pd.read_csv(args.inDesign)
    
    # Convert to a tall dataset
    countsTallDF = pd.melt(countsWideDF, id_vars=['featureID'], var_name='sampleID', value_name='APN')
    
    # Send tall counts dataset and design file to SQL and merge
    countsTallDF.to_sql("countsTall", con, if_exists="replace")
    designDF.to_sql("designInfo", con, if_exists="replace")
    cur.execute("CREATE TABLE countsKey AS SELECT in1.*, in2.featureID, in2.APN "
                "FROM designInfo in1 INNER JOIN countsTall in2 "
                "ON in1.sampleID = in2.sampleID ;")
    
    # Flag by event and sample, based on user-defined (or default) parameters
    ## Get user-defined thresholds. If they don't exist then set to 0
    if args.userAPN :
        minAPN = float(args.userAPN)
    else :
        minAPN = 0.0
    if args.userProp :
        minProp = float(args.userProp)
    else :
        minProp = 0.5
        
    # Get user-defined group
    userGroup = args.userGroup
        
    # First, flag each individual observation per sample if above the user-defined APN threshold
    cur.execute("ALTER TABLE countsKey ADD flag_on INT;")
    if minAPN == 0:
        cur.execute("UPDATE countsKey SET flag_on = (CASE WHEN APN > 0 THEN '1' ELSE '0' END);")
    else:
        cur.execute("UPDATE countsKey SET flag_on = (CASE WHEN APN > ? THEN '1' ELSE '0' END);", (minAPN,))

    # Summarize detection per group
    cur.execute("CREATE TABLE countsFlagged AS SELECT featureID, "+userGroup+", avg(flag_on) AS perc_trt FROM countsKey GROUP BY featureID, "+userGroup+" ;")
#    countsFlagged = pd.read_sql("SELECT * FROM countsFlagged;",con)
#    countsKey = pd.read_sql("SELECT * FROM countsKey;",con)

    # Create a list of groups to iterate over
    groupList = designDF[userGroup].drop_duplicates(keep='first').tolist()
    counter = 1
    for group in range(0,len(groupList)) :
        groupName = groupList[group]
        cur.execute("CREATE TABLE flagTemp AS SELECT featureID, perc_trt FROM countsFlagged "
                    "WHERE " + userGroup + " = ? ORDER BY featureID;", (groupName, ))
        cur.execute("ALTER TABLE flagTemp ADD flag_" + groupName + "_on"+str(int(minAPN))+" INT ;")
        if minProp == 0:  # min Prop = 0
            cur.execute("UPDATE flagTemp SET flag_" + groupName + "_on"+str(int(minAPN))+"  = (CASE WHEN perc_trt > 0 THEN '1' ELSE '0' END) ;")
        else:
            cur.execute("UPDATE flagTemp SET flag_" + groupName +
                            "_on"+str(int(minAPN))+" = (CASE WHEN perc_trt > ? THEN '1' ELSE '0' END) ;", (minProp, ))
        if counter == 1: # if first group value, then we create the output detection flags summary with group 1
            cur.execute("CREATE TABLE flagDtct AS SELECT featureID, flag_" + groupName + "_on"+str(int(minAPN))+" FROM flagTemp ; ")
            cur.execute("DROP TABLE flagTemp;")
            counter = counter + 1
        else:   # if not group 1, merge flagTemp with flags summary
            cur.execute("CREATE TABLE flagDtct2 AS SELECT in1.*, in2.flag_"
                        + groupName + "_on"+str(int(minAPN))+" FROM flagDtct in1 INNER JOIN flagTemp in2 "
                                      "ON in1.featureID = in2.featureID "
                                      "ORDER BY in1.featureID ;")
            cur.execute("DROP TABLE flagDtct;")
#            flagTemp = pd.read_sql("SELECT * FROM flagTemp;",con)
            cur.execute("DROP TABLE flagTemp;")
            cur.execute("ALTER TABLE flagDtct2 RENAME TO flagDtct ;")


    # Convert flags to tall dataset
    flagTall = pd.melt(pd.read_sql("SELECT * FROM flagDtct;", con), id_vars=['featureID'], var_name=userGroup, value_name='flag')
    flagTall[userGroup] = flagTall[userGroup].str.strip("flag_").str.strip("_on"+str(int(minAPN)))

    # Send tall counts dataset and design file to SQL and merge
    flagTall.to_sql("flagTall", con, if_exists="replace")
    cur.execute("CREATE TABLE flagKey AS SELECT in1.*, in2.featureID, in2.flag "
                "FROM designInfo in1 INNER JOIN flagTall in2 "
                "ON in1."+userGroup+" = in2."+userGroup+" ;")
#    flagKey = pd.read_sql("SELECT * FROM flagKey;",con)
    
    # Create other flags from sampleType flags
    if args.userFlags is not None:
        for other in args.userFlags:
            # Summarize detection per group
            cur.execute("CREATE TABLE otherFlagged AS SELECT featureID, "+other+", avg(flag) AS perc_trt FROM flagKey GROUP BY featureID, "+other+" ;")
#            otherFlagged = pd.read_sql("SELECT * FROM otherFlagged;",con)

            # Create a list of groups to iterate over
            otherList = designDF[other].drop_duplicates(keep='first').tolist()
            counter = 1
            for group in range(0,len(otherList)):
                groupName = otherList[group]
                cur.execute("CREATE TABLE otherTemp AS SELECT featureID, perc_trt FROM otherFlagged "
                            "WHERE " + other + " = ? ORDER BY featureID;", (groupName, ))
                cur.execute("ALTER TABLE otherTemp ADD flag_" + groupName + "_on"+str(int(minAPN))+" INT ;")
                if minProp == 0:  # min Prop = 0
                    cur.execute("UPDATE otherTemp SET flag_" + groupName + "_on"+str(int(minAPN))+"  = (CASE WHEN perc_trt > 0 THEN '1' ELSE '0' END) ;")
                else:
                    cur.execute("UPDATE otherTemp SET flag_" + groupName +
                                    "_on"+str(int(minAPN))+" = (CASE WHEN perc_trt > ? THEN '1' ELSE '0' END) ;", (minProp, ))
                if counter == 1: # if first group value, then we create the output detection flags summary with group 1
                    cur.execute("CREATE TABLE otherDtct AS SELECT featureID, flag_" + groupName + "_on"+str(int(minAPN))+" FROM otherTemp ; ")
                    cur.execute("DROP TABLE otherTemp;")
                    counter = counter + 1
                else:   # if not group 1, merge flagTemp with flags summary
                    cur.execute("CREATE TABLE otherDtct2 AS SELECT in1.*, in2.flag_"
                                + groupName + "_on"+str(int(minAPN))+" FROM otherDtct in1 INNER JOIN otherTemp in2 "
                                              "ON in1.featureID = in2.featureID "
                                              "ORDER BY in1.featureID ;")
                    cur.execute("DROP TABLE otherDtct;")
                    cur.execute("DROP TABLE otherTemp;")
                    cur.execute("ALTER TABLE otherDtct2 RENAME TO otherDtct ;")
            # Merge otherDtct with full table
            if args.userFlags[0] == other:
                cur.execute("CREATE TABLE allDtct AS SELECT in1.*, in2.* "
                            "FROM flagDtct in1 INNER JOIN otherDtct in2 "
                            "ON in1.featureID = in2.featureID ;")
            else:
                cur.execute("CREATE TABLE allDtct2 AS SELECT in1.*, in2.* "
                            "FROM allDtct in1 INNER JOIN otherDtct in2 "
                            "ON in1.featureID = in2.featureID ;")
                cur.execute("DROP TABLE allDtct;")
                cur.execute("ALTER TABLE allDtct2 RENAME TO allDtct;")
            
        detectedAllDF = pd.read_sql("SELECT * FROM allDtct;",con)

    
    # Create all sampleID on and off flags
    if args.userAll is not None:
        for allsID in args.userAll:
            # Summarize detection per group
            cur.execute("CREATE TABLE allsFlagged as SELECT featureID, "+allsID+", avg(flag_on) AS perc_trt FROM countsKey GROUP BY featureID, "+allsID+" ;")
#            allFlagged = pd.read_sql("SELECT * FROM allsFlagged;",con)
        
            # Create a list of groups to iterate over
            allList = designDF[allsID].drop_duplicates(keep='first').tolist()
            counter = 1

            for group in range(0,len(allList)):
                groupName = allList[group]
                cur.execute("CREATE TABLE allsTemp AS SELECT featureID, perc_trt FROM allsFlagged "
                            "WHERE " + allsID + " = ? ORDER BY featureID;", (groupName, ))
                cur.execute("ALTER TABLE allsTemp ADD flag_" + groupName + "_allsID_on"+str(int(minAPN))+" INT ;")
                cur.execute("ALTER TABLE allsTemp ADD flag_" + groupName + "_allsID_off"+str(int(minAPN))+" INT ;")
                # All sampleID's must be detected or not detected (perc_trt = 1 or 0)
                cur.execute("UPDATE allsTemp SET flag_" + groupName + "_allsID_on"+str(int(minAPN))+"  = (CASE WHEN perc_trt = 1 THEN '1' ELSE '0' END) ;")
                cur.execute("UPDATE allsTemp SET flag_" + groupName + "_allsID_off"+str(int(minAPN))+"  = (CASE WHEN perc_trt = 0 THEN '1' ELSE '0' END) ;")
                if counter == 1: # if first group value, then we create the output detection flags summary with group 1
                    cur.execute("CREATE TABLE allsDtct AS SELECT featureID, flag_" + groupName + "_allsID_on"+str(int(minAPN))+", flag_" + groupName + "_allsID_off"+str(int(minAPN))+" FROM allsTemp ; ")
                    cur.execute("DROP TABLE allsTemp;")
                    counter = counter + 1
                else:   # if not group 1, merge flagTemp with flags summary
                    cur.execute("CREATE TABLE allsDtct2 AS SELECT in1.*, in2.flag_"+groupName+
                                "_allsID_on"+str(int(minAPN))+", in2.flag_"+groupName +
                                "_allsID_off"+str(int(minAPN))+" FROM allsDtct in1 INNER JOIN allsTemp in2 "
                                              "ON in1.featureID = in2.featureID "
                                              "ORDER BY in1.featureID ;")
                    cur.execute("DROP TABLE allsDtct;")
                    cur.execute("DROP TABLE allsTemp;")
                    cur.execute("ALTER TABLE allsDtct2 RENAME TO allsDtct ;")
            # Merge otherDtct with full table
            cur.execute("CREATE TABLE allDtct2 AS SELECT in1.*, in2.* "
                        "FROM allDtct in1 INNER JOIN allsDtct in2 "
                        "ON in1.featureID = in2.featureID ;")
            cur.execute("DROP TABLE allDtct;")
            cur.execute("ALTER TABLE allDtct2 RENAME TO allDtct;")
        
        detectedAllDF = pd.read_sql("SELECT * FROM allDtct;",con)
        cols = [c for c in detectedAllDF.columns if c[:10] != "featureID:"]
        detectedAllDF = detectedAllDF[cols]
    
    # Write output flags file
    detectFlagDF = pd.read_sql("SELECT * FROM flagDtct;", con)
    if args.userFlags is None and args.userAll is None:
        detectFlagDF.to_csv(args.outFlags, encoding='utf-8', index=False)
    else:
        detectedAllDF.to_csv(args.outFlags, encoding='utf-8', index=False)

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

