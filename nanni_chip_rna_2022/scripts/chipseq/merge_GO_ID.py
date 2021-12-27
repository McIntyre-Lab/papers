#!/usr/bin/env python

import argparse
import pandas as pd
import sqlite3


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge in Gene Ontology IDs based on gene IDs")

    # Input data
    parser.add_argument("-f", "--flags", dest="inFlags", required=True, help="Input CSV file with NO HEADER of reference gene_id in first column and matching PB gene_id in the second coloumn")
    parser.add_argument("-g", "--go", dest="inGO", required=True, help="Input CSV file with Gene to GO IDs (columns: biological_process, molecular_function, cellular_component)")
    parser.add_argument("--name-flag", dest="inNF", required=True, help="Name of gene ID column in flag file to merge on")
    parser.add_argument("--name-go", dest="inNG", required=True, help="Name of gene ID column in GO file to merge on")

    # Output data
    parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Output prefix for flags with go terms")

    args = parser.parse_args()
    return args

def main():
    
    # Get list of flags
    flagDF = pd.read_csv(args.inFlags)
    flagName = args.inNF
    
    # Get GO ID file
    goDF = pd.read_csv(args.inGO)
    goName = args.inNG
    
    # Merge GO terms by gene_id
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    flagDF.to_sql("flag", con, if_exists="replace")
    goDF.to_sql("go", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.biological_process, in2.molecular_function, in2.cellular_component "
                "FROM flag in1 LEFT JOIN go in2 "
                "ON in1."+flagName+" = in2."+goName+";")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])

    mergeDF.to_csv(args.outPrefix+".goID.csv", index=False)
    
    withGO = mergeDF[(~mergeDF['biological_process'].isna()) | (~mergeDF['molecular_function'].isna()) | (~mergeDF['cellular_component'].isna())]
    noGO = mergeDF[(mergeDF['biological_process'].isna()) & (mergeDF['molecular_function'].isna()) & (mergeDF['cellular_component'].isna())]
    print("Genes without GO ID = " + str(len(noGO)))
    
    withGO.to_csv(args.outPrefix+".goID.withGO.csv",index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

