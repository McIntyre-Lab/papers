#!/usr/bin/env python
# Imports a BED file for extracting transcript FASTA sequence and converts it to a BED file of exon sequences

# Built-in packages
import argparse
from operator import itemgetter

# Add-on packages
import pandas as pd
import numpy as np
import csv

# Arguments
parser = argparse.ArgumentParser(description="Converts a 12-column BED file of transcripts to a 6-column BED file"
                                             "of exons")
parser.add_argument("--input", dest="inputBED", action='store', required=True,
                    help="Input 12-column BED file. Columns 11 and 12 should refer to block (exon) sizes and their"
                         "genomic position relative to the genomic start of the transcript")
parser.add_argument("--output", dest="outputBED", action='store', required=True,
                    help="Output 6-column BED file of exons for extracting FASTA sequence")

args = parser.parse_args()

with open(args.outputBED, 'w') as outputFile:
    xsBED = pd.read_csv(args.inputBED, sep="\t", header=None,
                        names=['chr','start','stop','transcriptID','score','strand','thickStart','thickStop',
                               'itemRGB','blockCount','blockSize','blockStart'])
    # Set up lists for creating exon BED file
    exonChrList = []
    exonStartList = []
    exonStopList = []
    exonIDList = []
    exonScoreList = []
    exonStrandList = []
    for index, row in xsBED.iterrows():
        if np.isnan(row['start']) == False and row['blockCount'] > 0 :
            exSizes=[]
            exStarts = []
            exCounter = 1
            numExons = int(row['blockCount'])
            # Put block starts and sizes into lists to iterate over
            if "," in row['blockSize'] :
                exSizes=row['blockSize'].split(",")
                exStarts = row['blockStart'].split(",")
            else:
                exSizes.append(row['blockSize'])
                exStarts.append(row['blockStart'])
            # Iterate over block lists and append exon info to lists
            for ex in range(0,numExons):
                exonStart=int(row['start']) + int(exStarts[ex])
                exonStop=exonStart + int(exSizes[ex])
                exonID=str(row['transcriptID']) + ":" + str(exCounter)
                # Append to exon Lists
                exonChrList.append(row['chr'])
                exonStartList.append(exonStart)
                exonStopList.append(exonStop)
                exonScoreList.append(".")
                exonStrandList.append(row['strand'])
                exonIDList.append(exonID)
                exCounter=exCounter+1
    # Make a new dataframe from exon lists and write dataframe to output
    exonBED = pd.DataFrame({'chr':exonChrList, 'start':exonStartList, 'stop': exonStopList,'exon_id': exonIDList,
                            'score': exonScoreList,'strand': exonStrandList })
    exonBED.to_csv(outputFile, encoding='utf-8', index=False, sep='\t', header=None, quoting=csv.QUOTE_NONE)
