#!/usr/bin/env python

import argparse
import pandas as pd
import sqlite3

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge transcript-level junction, fragment, and exon region variables")

    # Input data
    parser.add_argument("-j", "--junctions", dest="inJunc", required=True, help="CSV file of transcript-level junction variables")
    parser.add_argument("-f", "--fragments", dest="inFrag", required=True, help="CSV file of transcript-level fragment variables")
    parser.add_argument("-e", "--exon-regions", dest="inER", required=True, help="CSV file of transcript-level exon region variables")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file of merged transcript-level junction, fragment, and exon region variables")

    args = parser.parse_args()
    return args

def main():
    # Get input file of transcript-level junction variables
    collapseJuncDF = pd.read_csv(args.inJunc, low_memory=False)

    # Get input file of transcript-level fragment variables
    collapseFragDF = pd.read_csv(args.inFrag, low_memory=False) 
    
    # Get input file of transcript-level exon region variables
    collapseExonRegionDF = pd.read_csv(args.inER, low_memory=False) 
    
    # Verify that all input dataframes are the same length
    if (len(collapseJuncDF) != len(collapseFragDF)) | (len(collapseJuncDF) != len(collapseExonRegionDF)):
        print("WARNING: Input transcript-level files are different lengths")

    # Merge junction variables with fragment variables
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    collapseJuncDF .to_sql("junction", con, if_exists="replace")
    collapseFragDF.to_sql("fragment", con, if_exists="replace")
    cur.execute("CREATE TABLE mergeJuncFrag AS SELECT in1.gene_id, in1.transcript_id, "
                "in1.junctionID_order, in1.flag_monoexon_transcript, in2.fragmentID_order, "
                "in2.num_fragment, in2.chr, in2.start, in2.end, in2.transcript_length, in2.fragmentLength_order "
                "FROM junction in1 LEFT JOIN fragment in2 "
                "ON in1.transcript_id = in2.transcript_id")
    del(collapseJuncDF, collapseFragDF)
    # Merge exon region variables with previous merge of junction and fragment variables
    collapseExonRegionDF.to_sql("exonRegion", con, if_exists="replace")
    cur.execute("CREATE TABLE mergeAll AS SELECT in1.*, in2.exonRegionID_order, "
                "in2.exonRegionLength_order, in2.num_exonRegion "
                "FROM mergeJuncFrag in1 LEFT JOIN exonRegion in2 "
                "ON in1.transcript_id = in2.transcript_id")
    mergeAllDF = pd.read_sql("SELECT * FROM mergeAll", con)
    con.close()
    del(collapseExonRegionDF)
    
    # Output CSV of transcripts with splice group information (NOT collapsed yet)
    mergeAllDF.to_csv(args.outFile, index=False)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

