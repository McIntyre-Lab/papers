#!/usr/bin/env python

import argparse
import pandas as pd
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Subset Event Analysis (EA) annotation files using a transcript exclusion list")

    # Input data
    parser.add_argument("-i", "--input-annotation", dest="inAnnot", required=True, help="EA annotation CSV file (must contain columns: transcript_id and annotation_frequency)")
    parser.add_argument("-v", "--input-index-variable", dest="inVar", required=True, help="Index variable for input EA annotation CSV file (eg. fragment_id, fusion_id, event_id, etc.)")
    parser.add_argument("-l", "--exclusion-list", dest="inList", required=True, help="List of transcripts to be excluded from the input EA file (no header, 1 transcript_id per line)")

    # Output data
    parser.add_argument("-o", "--outFile", dest="outFile", required=True, help="Output file path and name for subset EA annotation CSV")

    args = parser.parse_args()
    return args

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
    # Split variable by some character like '|' or ',' and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def main():
    # Get input EA annotation file and exclusion list
    annotDF = pd.read_csv(args.inAnnot, low_memory=False)
    excludeDF = pd.read_csv(args.inList, names=['transcript_id'], low_memory=False)
    
    # Set index variable
    indexVar = args.inVar
    if indexVar not in annotDF.columns:
        print("ERROR: Index variable name provided is not in EA input CSV file")
        sys.exit()
    if len(annotDF[indexVar]) != annotDF[indexVar].nunique():
        # index is not unique - annotation file may already be split by transcript
        #     (therefore no "|" should be present in transcript_id column)
        if annotDF['transcript_id'].str.find("|").max() > -1: # No "|" found
            print("ERROR: Index variable provided is not unique in EA input CSV file")
            sys.exit()
        else:
            # Subset annotation file by excluding the transcripts in the list provided
            subsetDF = annotDF[~annotDF['transcript_id'].isin(excludeDF['transcript_id'])]
    else:
        # index variable is unique therefore transcripts are potentially concatenated by "|"
        # Split transcript_id by "|"
        splitDF = split_column_by_sep(df=annotDF,col_name='transcript_id',sep="|")
        # Subset the split annotation file by excluding the transcript in the list provided
        splitSubsetDF = splitDF[~splitDF['transcript_id'].isin(excludeDF['transcript_id'])]
        # Make an aggregation function dictionary where all columns select the 'first'
        #     except for the variable grouped by (indexVar) and transcript_id which will be
        #     concatenated with "|" to match the format of the input annotation file
        #     (Concatenated lists will no longer contain excluded transcript_id values)
        funcDict = {}
        for col in splitSubsetDF.columns:
            if col == indexVar:
                continue
            elif col == 'transcript_id':
                funcDict.update({col:[lambda x: "|".join(x.map(str))]})
            else:
                funcDict.update({col:['first']})
        # Aggregate columns using above functions
        subsetDF = splitSubsetDF.groupby(indexVar).agg(funcDict).reset_index()
        # Fix column values
        subsetDF.columns = subsetDF.columns.droplevel(1)
        
    # Output subset annotation file
    subsetDF.to_csv(args.outFile,index=False)
    

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
