#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import sqlite3

# Import custon functions for FSM reduction
import FSM_consolidation_functions as FSM

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get transcript-level junction variables")

    # Input data
    parser.add_argument("-i", "--input-transcripts", dest="inXcrptGene", required=True, help="CSV file with unique pairs of transcript_id to gene_id")
    parser.add_argument("-j", "--junctions", dest="inJunc", required=True, help="CSV file of junction annotations (*_annotated_junctions.csv)")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file of transcript-level junction variables")

    args = parser.parse_args()
    return args

def main():
    # Get input file of transcript_id to gene_id
    #     (will be used to get gene_id for monoexon transcripts)
    xcrptGene = pd.read_csv(args.inXcrptGene, low_memory=False)
    
    # Get input file of junctions (will only include spliced transcripts)
    #      junction_id is not unique (pair of junction_id,transcript_id is unique)
    juncDF = pd.read_csv(args.inJunc, low_memory=False)
    # Drop empty junction_id values if present
    juncDF = juncDF.dropna()
    # Drop exon ID junction_id and set junction_coordinates as junction_id
    #     (matches the junction_id values in *_junctions_full_annotation.csv)
    juncDF = juncDF.drop(columns=['junction_id'])
    juncDF = juncDF.rename(columns={'junction_coordinates':'junction_id'})
    
    # Get donor_stop and acceptor_start from coordinates in junction_id
    juncDF['donor_stop'] = juncDF['junction_id'].str.split(":").str[1].astype(int)
    juncDF['acceptor_start'] = juncDF['junction_id'].str.split(":").str[2].astype(int)
    
    # Count junctions
    print("{} unique junctions present in annotation".format(juncDF['junction_id'].nunique()))
    
    # Make transcript-level junction variables (including junctionID_order)
    #      by collapsing by transcript_id
    # Sort junctions by transcript_id and donor_stop
    juncDF = juncDF.sort_values(['transcript_id','donor_stop'])
    collapseJuncDF = FSM.collapse_junctions(juncDF)
    del(juncDF)

    # Merge transcript-level junction (multiexon transcripts only) and
    #     gene_id values (all transcripts including multi-exon) variables to get
    #     mono-exon transcripts
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    xcrptGene.to_sql("gene", con, if_exists="replace")
    collapseJuncDF.to_sql("junction", con, if_exists="replace")
    # Merge transcript level junction variables and gene_id by transcript_id
    cur.execute("CREATE TABLE merge1 AS SELECT in1.gene_id, in1.transcript_id, in2.junctionID_order "
                "FROM gene in1 LEFT JOIN junction in2 "
                "ON in1.transcript_id = in2.transcript_id")
    collapseJuncXcrptDF = pd.read_sql("SELECT * FROM merge1", con)
    # Verify merge is good
    if len(collapseJuncXcrptDF) != len(xcrptGene):
        print("WARNING: Unexpected merge")
    if len(collapseJuncXcrptDF) != len(collapseJuncDF) + len(collapseJuncXcrptDF[collapseJuncXcrptDF['junctionID_order'].isna()]):
        print("WARNING: Unexpected merge")
    del(collapseJuncDF, xcrptGene)

    # Flag monoexon transcripts (no junctionID_order present)
    collapseJuncXcrptDF['flag_monoexon_transcript'] = np.where(collapseJuncXcrptDF['junctionID_order'].isna(),1,0)

    # Transcripts with no junctionID_order will be set to transcript_id value
    collapseJuncXcrptDF = collapseJuncXcrptDF.fillna(-1)
    collapseJuncXcrptDF.loc[collapseJuncXcrptDF['junctionID_order']==-1,'junctionID_order'] = collapseJuncXcrptDF['transcript_id']

    # Output transcript-level file
    collapseJuncXcrptDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

