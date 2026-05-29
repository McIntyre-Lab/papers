#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import csv
import sqlite3
#import sys


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Correct the gene_id values of the GTF output from SQANTI3 QC to match the gene_id values associated with the transcripts in the classification file")

    # Input data
    parser.add_argument("-c", "--classification", dest="inClass", required=True, help="Input SQANTI classification file")
    parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="Input GTF file to be modified")

    # Output data
    parser.add_argument("--output-classification", dest="outClass", required=False, help="Output TSV file of classification file lines without associated genes")
    parser.add_argument("--output-gtf", dest="outGTF", required=False, help="Output GTF file of isoforms that did not get associated genes assigned in classification file")
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output GTF file name for isoform exons with associated gene_id (file is UNSORTED)")

    args = parser.parse_args()
    return args

def main():

    # Get SQANTI QC classification output 
    classDF = pd.read_csv(args.inClass, sep="\t", low_memory=False)
    
    # Output classification file lines that did not get an associated gene if argument used
    if args.outClass is not None:
        noAssocGene = classDF[classDF['associated_gene'].isna()]
        noAssocGene.to_csv(args.outClass,sep="\t",index=False)

#    tempDF = classDF[classDF['associated_gene'].isna()]
    
    # Select isoform and associated_gene values
    classDF = classDF[['isoform','associated_gene']]
    
    # Get GTF file and count total lines
    gtf = pd.read_csv(args.inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], sep="\t",low_memory=False)
    print("Total lines in original GTF= "+str(len(gtf)))

    # Extract gene_id from attribute column (this is actually the isoform name)
    gtf['gene_id'] = gtf['attribute'].str.split(";").str[1].str.split("\"").str[1]
    
    # Merge associated_gene where gene_id = isoform
    con = sqlite3.connect(':memory:')
    cur = con.cursor()
    classDF.to_sql("classification", con, if_exists="replace")
    gtf.to_sql("gtf", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.associated_gene "
                "FROM gtf in1 LEFT JOIN classification in2 "
                "ON in1.gene_id = in2.isoform")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])
    con.close()
    
    # Check for proper merge
#    if len(mergeDF[mergeDF['associated_gene'].isna()]) > 0:
#        sys.exit("!!! ERROR : UNEXPECTED MERGE RESULTING IN MISSING associated_gene")
    
    # Output GTF lines that did not get an associated gene
    testGTF = mergeDF[mergeDF['associated_gene'].isna()]
    testGTF[['chr','source','feature','start','end','score','strand','frame',
             'attribute']].to_csv(args.outGTF,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)
    
    # Create new attribute column with associated_gene as gene_id
    mergeDF['new_attribute'] = mergeDF['attribute'].str.split(";").str[0] + "; gene_id \"" + mergeDF['associated_gene'] + "\";"

    # new_attribute will equal original attribute if there is no associated gene
    mergeDF['new_attribute'] = np.where(mergeDF['associated_gene'].isna(),mergeDF['attribute'],mergeDF['new_attribute'])

    # Output GTF file with associated_gene in new attribute
    mergeDF[['chr','source','feature','start','end','score','strand','frame',
             'new_attribute']].to_csv(args.outFile,sep="\t",index=False,header=False,
             doublequote=False,quoting=csv.QUOTE_NONE)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

