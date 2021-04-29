#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Count transcripts per gene")

    # Input arguments
    parser.add_argument("-g", "--gtf", dest="inGTF", required=True, help="GTF file where transcript_id and gene_id are the first two attributes in the attributes column (order does not matter)")
 
    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file of counts")

    args = parser.parse_args()
    return args


def main():
    # Get GTF file and count total lines
    gtf = pd.read_csv(args.inGTF,names=['chr','source','feature','start','end','score',
                                        'strand','frame','attribute'], dtype=str, sep="\t",low_memory=False)

     # Get attribute labels 1 and 2
    gtf['attribute_name_1'] = gtf['attribute'].str.split(" ").str[0]
    gtf['attribute_name_2'] = gtf['attribute'].str.split(" ").str[2]
    
    # Extract transcript_id and gene_id from attribute column
    gtf['transcript_id'] = np.where(gtf['attribute_name_1']=="transcript_id",
                                       gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                       np.where(gtf['attribute_name_2']=="transcript_id",
                                                gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))
    gtf['gene_id'] = np.where(gtf['attribute_name_1']=="gene_id",
                                       gtf['attribute'].str.split(";").str[0].str.split("\"").str[1],
                                       np.where(gtf['attribute_name_2']=="gene_id",
                                                gtf['attribute'].str.split(";").str[1].str.split("\"").str[1],"oops"))

    # Open output file
    outFile = open(args.outFile,'w')

    # Get transcripts per gene
    outFile.write("\n{} unique transcripts in {} unique genes".format(gtf['transcript_id'].nunique(),gtf['gene_id'].nunique()))
    xcrptPerGene = gtf.groupby('gene_id')['transcript_id'].nunique().reset_index().rename(columns={'transcript_id':'transcript_per_gene'})
    freqXcrptPerGene = xcrptPerGene['transcript_per_gene'].value_counts()
    outFile.write("\nBins of transcripts per gene frequencies:\ntranscript_per_gene_bin\tfrequency")
    outFile.write("\n\t\t1\t\t\t\t\t{}".format(freqXcrptPerGene[1]))
    outFile.write("\n\t\t2\t\t\t\t\t{}".format(freqXcrptPerGene[2]))
    outFile.write("\n\t\t3\t\t\t\t\t{}".format(freqXcrptPerGene[3]))
    outFile.write("\n\t\t[4-6)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(4,6) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[6-8)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(6,8) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[8-10)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(8,10) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[10-12)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(10,12) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[12-14)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(12,14) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[14-16)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(14,16) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[16-18)\t\t\t\t{}".format(freqXcrptPerGene[[c for c in range(16,18) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\n\t\t[18-{}]\t\t\t\t{}".format(xcrptPerGene['transcript_per_gene'].max(),freqXcrptPerGene[[c for c in range(18,xcrptPerGene['transcript_per_gene'].max()+1) if c in freqXcrptPerGene.index]].sum()))
    outFile.write("\nAll transcript per gene value frequencies:\n{}".format(freqXcrptPerGene.sort_index().to_string()))

    outFile.close()
    
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

