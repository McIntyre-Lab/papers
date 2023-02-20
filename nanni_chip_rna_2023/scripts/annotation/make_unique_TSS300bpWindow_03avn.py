#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sqlite3
import gffutils
import argparse
import sys
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Obtain coordinates of all annotated 5'UTR, 3'UTR, and TSS (1kb region centered on 5'UTR start) from reference GTF")

    # Input data
    parser.add_argument("-i", "--inGTF", dest="inGTF", required=True, help="Input GTF, if there is no corresponding database (*.db) file it will be created")
    parser.add_argument("-g", "--genome", dest="inChrom", required=True, help="Genome chromosome sizes, format: chr\tsize")

    # Output data
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temporary SQL database file to reduce memory requiredfor SQL merge", default=":memory:")
    parser.add_argument("-o", "--output-prefix", dest="outPrefix", required=True, help="Output prefix for the bed file of unique values and CSV for overlapID to featureIDs")

    args = parser.parse_args()
    return args

def main():

    # Make chromosome sizes DF
    chromDF = pd.read_csv(args.inChrom, sep="\t", names=['chr','size'])
    
    # Initialize TSS dataframe
    tssDF = pd.DataFrame(columns=['chrom','start','end','FBgn','FBtr','featureType'])

    # Open GTF
    if not os.path.isfile(args.inGTF+".db"):
        db =  gffutils.create_db(args.inGTF, dbfn=args.inGTF+".db", force=True, force_gff=False)
    else:
        db = gffutils.FeatureDB(args.inGTF+".db")
        
    # Use transcript starts to set TSS1kbWindow (+/- 500bp)
    for i in db.features_of_type('transcript'):
        gene = str(i.attributes['gene_id']).strip('[\'').strip('\']')
        trans = str(i.attributes['transcript_id']).strip('[\'').strip('\']')
        if i.strand == "+":
            newStart = i.start - 151
            newEnd = i.start + 150
            if newStart < 0 : newStart = 0
            if newEnd > chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] :
                newEnd = chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0]
        else :
            newStart = i.end - 151
            newEnd = i.end + 150
            if newStart < 0 : newStart = 0
            if newEnd > chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] :
                newEnd = chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] 
        tssDF.loc[len(tssDF)]=[str(i.chrom)]+[str(newStart)]+[str(newEnd)]+[gene]+[trans]+["TSS"]
        
    # Get featureID names
    newDF = tssDF[['chrom','start','end']].drop_duplicates()
    newDF['featureID'] = "TSS300bpWindow_" + newDF['chrom'].astype(str) + "_" + newDF['start'].astype(str) + "_" + newDF['end'].astype(str)

    # Merge names with all features
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    
    tssDF.to_sql("features", con, if_exists="replace")
    newDF.to_sql("names", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.featureID "
                "FROM features in1 INNER JOIN names in2 "
                "ON in1.chrom = in2.chrom AND in1.start = in2.start AND in1.end = in2.end ")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])
    
    # Check length of mergeDF - should equal the length of the input
    if len(mergeDF) != len(tssDF):
        print("ERROR : Improper merge of unique values with input")
        sys.exit()

    # Output bed file of unique features
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
