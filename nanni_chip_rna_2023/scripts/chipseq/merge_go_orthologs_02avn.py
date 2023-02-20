#!/usr/bin/env python

import argparse
import pandas as pd
import sqlite3

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Merge in Gene Ontology IDs based on gene IDs")

    # Input data
    parser.add_argument("-g", "--go", dest="inGO", required=True, help="Input CSV file with Gene to GO IDs (columns: biological_process, molecular_function, cellular_component)")
    parser.add_argument("-o", "--ortho", dest="inOrtho", required=True, help="Input CSV file of reformated orthologous genes")

    # Output data
    parser.add_argument("--output", dest="outFile", required=True, help="Output CSV file with melanogaster and simulans GO annotations based on orthologs")
    
    args = parser.parse_args()
    return args

def main():

    # Get input files
    goDF = pd.read_csv(args.inGO)
    orthoDF = pd.read_csv(args.inOrtho, names=['mel_geneID','mel_geneSymbol','mel_chrom','mel_start','mel_end','mel_strand','sim_geneID','sim_geneSymbol','sim_chrom','sim_start','sim_end','sim_strand','orthoDB_groupID'])
    
    print("Number of mel genes with GO annotations = "+str(len(goDF)))
    
    # Merge GO terms by gene_id
    con = sqlite3.connect(":memory:")
    cur = con.cursor()

    orthoDF.to_sql("ortho", con, if_exists="replace")
    goDF.to_sql("go", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in2.sim_geneID, in2.sim_geneSymbol, in2.mel_geneID, in2.mel_geneSymbol, in1.* "
                "FROM go in1 LEFT JOIN ortho in2 "
                "ON in1.FBgn = in2.mel_geneID ;")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])

    mergeDF = mergeDF[~mergeDF['sim_geneID'].isna()].drop(columns=['FBgn','symbol'])
    
    # Count each gene and drop all paralogs
    mergeDF['mel_cnt'] = mergeDF.groupby('mel_geneID')['mel_geneID'].transform('count')
    mergeDF['sim_cnt'] = mergeDF.groupby('sim_geneID')['sim_geneID'].transform('count')
    allParDF = mergeDF[(mergeDF['mel_cnt']>1) & (mergeDF['sim_cnt']>1)].drop(columns=['mel_cnt','sim_cnt'])
    noParDF = mergeDF[(mergeDF['mel_cnt']==1) & (mergeDF['sim_cnt']==1)].drop(columns=['mel_cnt','sim_cnt'])
    print("Number of paralogous pairs with GO annnotations = "+str(len(allParDF)))
    print("Number of unique orthologous pairs with GO annotations = "+str(len(noParDF)))

    noParDF.to_csv(args.outFile, index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()


