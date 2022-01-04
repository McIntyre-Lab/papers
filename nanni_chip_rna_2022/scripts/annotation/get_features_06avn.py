#!/usr/bin/env python3

import pandas as pd
import gffutils
import argparse
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Obtain coordinates of all annotated 5'UTR, 3'UTR, and TSS (1kb region centered on 5'UTR start) from reference GTF")

    # Input data
    parser.add_argument("-i", "--inGTF", dest="inGTF", required=True, help="Input GTF, if there is no corresponding database (*.db) file it will be created")
    parser.add_argument("-g", "--genome", dest="inChrom", required=True, help="Genome chromosome sizes, format: chr\tsize")

    # Output data
    parser.add_argument("-o", "--output-dir", dest="outDir", required=True, help="Output directory for unsorted bed files of each feature type")

    args = parser.parse_args()
    return args

def main():

    # Make chromosome sizes DF
    chromDF = pd.read_csv(args.inChrom, sep="\t", names=['chr','size'])
    
    # Open output files
    outTSS = open(args.outDir+'/all_TSS1kbWindow.bed','w')
    out5UTR = open(args.outDir+'/all_5UTR.bed','w')
    out3UTR = open(args.outDir+'/all_3UTR.bed','w')
    
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
            newStart = i.start - 501
            newEnd = i.start + 500
            if newStart < 0 : newStart = 0
            if newEnd > chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] :
                newEnd = chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0]
        else :
            newStart = i.end - 501
            newEnd = i.end + 500
            if newStart < 0 : newStart = 0
            if newEnd > chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] :
                newEnd = chromDF[chromDF['chr']==str(i.chrom)]['size'].values[0] 
        outTSS.write(str(i.chrom) + "\t" + str(newStart) + "\t" + str(newEnd) + "\t" + gene + "\t" + trans + "\tTSS\n")

    # Get all 5'UTR entries - Note: UTR can be split across multiple exons
    for i in db.features_of_type('5UTR') :
        gene = str(i.attributes['gene_id']).strip('[\'').strip('\']')
        trans = str(i.attributes['transcript_id']).strip('[\'').strip('\']')
        out5UTR.write(str(i.chrom) + "\t" + str(i.start-1) + "\t" + str(i.stop) + "\t" + gene + "\t" + trans + "\t5UTR\n" )

    ## Get all 3'UTR entries - Note: UTR can be split across multiple exons
    for i in db.features_of_type('3UTR') :
        gene = str(i.attributes['gene_id']).strip('[\'').strip('\']')
        trans = str(i.attributes['transcript_id']).strip('[\'').strip('\']')
        out3UTR.write(str(i.chrom) + "\t" + str(i.start-1) + "\t" + str(i.stop) + "\t" + gene + "\t" + trans + "\t3UTR\n")        
    outTSS.close()    
    out5UTR.close()
    out3UTR.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

