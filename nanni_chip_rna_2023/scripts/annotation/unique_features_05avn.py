#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sqlite3
import argparse
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in a bed file and outputs a bed of uniq values as well as a CSV file connecting overlapID to featureIDs")

    # Input data
    parser.add_argument("-b", "--bedFile", dest="inBed", required=True, help="Bed file with columns: chr start end FBgn FBtr featureType")
    
    # Output data
    parser.add_argument("-o", "--outPrefix", dest="outPrefix", required=True, help="Output prefix for the bed file of unique values and CSV for overlapID to featureIDs")
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temporary SQL database file to reduce memory requiredfor SQL merge", default=":memory:")
    
    args = parser.parse_args()
    return args

def main():

    # Import bed file of features
    inDF = pd.read_csv(args.inBed, sep="\t", names=['chrom','start','end','FBgn','FBtr','featureType'])
    
    # Check all featureIDs are of the same type
    if len(inDF['featureType'].unique()) > 1 :
        print("ERROR : Not all features of the same type")
        sys.exit()
    else :
        featType = inDF['featureType'].unique()[0]
    
    # Get featureID names
    newDF = inDF[['chrom','start','end']].drop_duplicates()
    newDF['featureID'] = featType + "_" + newDF['chrom'].astype(str) + "_" + newDF['start'].astype(str) + "_" + newDF['end'].astype(str)
    
    # Merge names with all features
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    
    inDF.to_sql("features", con, if_exists="replace")
    newDF.to_sql("names", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.featureID "
                "FROM features in1 INNER JOIN names in2 "
                "ON in1.chrom = in2.chrom AND in1.start = in2.start AND in1.end = in2.end ")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])

    # Check length of mergeDF - should equal the length of the input
    if len(mergeDF) != len(inDF):
        print("ERROR : Improper merge of unique values with input")
        sys.exit()

    newDF[['chrom','start','end','featureID']].to_csv(args.outPrefix+"_unique.bed", sep="\t", index=False, header=False)

    # Count genes and transcripts per featureID and flag multigene
    mergeDF['num_FBgn'] = mergeDF.groupby('featureID')['FBgn'].transform('nunique')
    mergeDF['flag_multigene'] = np.where(mergeDF['num_FBgn']>1,"1","0")
    mergeDF['num_FBtr'] = mergeDF.groupby('featureID')['FBtr'].transform('nunique')
    
    # Add annotation frequency information (similar to Event Analysis)
    # Must make sure multigene transcript/gene sums are correct
    mergeDF['num_xcrpt_per_gene'] = mergeDF.groupby('FBgn')['FBtr'].transform('nunique')
    tempSum = mergeDF[['featureID','FBgn','num_xcrpt_per_gene']].drop_duplicates().groupby(['featureID','FBgn'])['num_xcrpt_per_gene'].sum().reset_index()
    tempSum['num_xcrpt_per_all_genes'] = tempSum.groupby('featureID')['num_xcrpt_per_gene'].transform('sum')
    tempSum.to_sql("tempSums", con, if_exists="replace")
    mergeDF.to_sql("newMerge", con, if_exists="replace")
    cur.execute("CREATE TABLE key AS SELECT in1.*, in2.num_xcrpt_per_all_genes "
                "FROM newMerge in1 LEFT JOIN tempSums in2 "
                "ON in1.featureID = in2.featureID AND in1.FBgn = in2.FBgn ")
    sumDF = pd.read_sql("SELECT * FROM key", con).drop(columns=['index'])
    sumDF['annotation_frequency'] = np.where(sumDF['num_FBtr']==1,"Unique",
           np.where(sumDF['num_FBtr'] < sumDF['num_xcrpt_per_all_genes'],"Common",
                    np.where(sumDF['num_FBtr'] == sumDF['num_xcrpt_per_all_genes'],"Constitutive","-1")))
    sumDF = sumDF.drop(columns=['num_xcrpt_per_gene','num_xcrpt_per_all_genes'])

    # Output annotation file
    sumDF[['num_FBtr','num_FBgn']] = sumDF[['num_FBtr','num_FBgn']].astype('str')
    mergeDF.to_csv(args.outPrefix+"_xcrpt_annotation.csv", index=False)

    # Get unique featureID annotation (set multigene/multitranscript FBgn and FBtr values to NaN)
    uniqDF = sumDF.groupby('featureID').first().reset_index()
    uniqDF.loc[uniqDF['flag_multigene']=="1",'FBgn'] = np.nan
    uniqDF.loc[uniqDF['num_FBtr']!="1",'FBtr'] = np.nan
    uniqDF[['featureID','chrom','start','end','FBgn','FBtr','featureType','num_FBgn','flag_multigene','num_FBtr','annotation_frequency']].to_csv(args.outPrefix+"_unique.csv", index=False, na_rep='')

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
