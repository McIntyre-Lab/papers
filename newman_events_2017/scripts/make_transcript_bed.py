#!/usr/bin/env python

#       DESCRIPTION: Export BED file of transcripts by exon

# Built-in packages
import argparse
import os

# Add-on packages
import gffutils
import itertools


def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Import database and give exported file name")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF or GTF file")
    parser.add_argument("--output", dest="outBED", action='store', required=True, help="Output BED file name")
    
    args = parser.parse_args()
    return args

def main():
    # Get the database
    db_fn=args.gffInput + '.db'
    db=gffutils.FeatureDB(db_fn)
    
    
    # Open output files for writing
    bedOut = open(args.outBED, 'w')

    xscripts = db.features_of_type('mRNA')
    for xs in xscripts:
        chrom=xs.chrom
        totalStart=xs.start-1
        totalStop=xs.stop
        strand=xs.strand
        name=xs.id
        startList = []
        lengthList = []
        exonCount=0
        exons = list(db.children(xs, featuretype='exon', order_by='start'))
        for ex in exons:
            exonStart=(ex.start-1)-totalStart
            exonLength=ex.end-(ex.start-1)
            startList.append(str(exonStart))
            lengthList.append(str(exonLength))
            exonCount=exonCount + 1
        #print exonCount
        #print startList
        #print lengthList
        startString=','.join(startList)
        lengthString=','.join(lengthList)
        xsEntry=[chrom, str(totalStart), str(totalStop), name, '.', strand, str(totalStart), str(totalStop), '255,0,0', str(exonCount), lengthString, startString]
        bedEntry = '\t'.join(xsEntry) + "\n"
        bedOut.write(bedEntry)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

