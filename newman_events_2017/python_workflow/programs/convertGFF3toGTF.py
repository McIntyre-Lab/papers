#!/usr/bin/env python
# Import GFF3 file and convert to an Ensembl-style GTF file for refSQANTI.py script

# Built-in packages
import argparse
from operator import itemgetter

# Add-on packages
#import gffutils
import pandas as pd
import csv

# Arguments
parser = argparse.ArgumentParser(description="Import database and give exported file name")
parser.add_argument("--input", dest="gtfInput", action='store', required=True,
                    help="Input GFF3 file, sorted by chromosome, start and stop positions")
parser.add_argument("--output", dest="outputFile", action='store', required=True, help="Output GTF file name")

args = parser.parse_args()

with open(args.outputFile, 'w') as outputFile:
    gffDF = pd.read_csv(args.gtfInput, sep="\t", header=None,
                        names=['chr','source','feature_type','start','stop','score','strand','frame','attributes'],
                        dtype=str)
    # Set up all my needed columns for the final dataFrame
    parentList = []
    idList = []
    chrList = []
    sourceList = []
    typeList = []
    startList = []
    stopList = []
    scoreList = []
    strandList = []
    frameList = []
    attributeList = []
    for index, row in gffDF.iterrows():
        # Ensembl GTF files need the following: gene, transcript, exon, CDS
        if row['feature_type'] == "gene": # If entry is a gene
            chrList.append(row['chr'])
            sourceList.append(row['source'])
            typeList.append(row['feature_type'])
            startList.append(row['start'])
            stopList.append(row['stop'])
            scoreList.append(row['score'])
            strandList.append(row['strand'])
            frameList.append(row['frame'])
            # Reformat attributes -- from GFF3, there is geneID and gene name only
            geneID=(row['attributes'].split(";")[0]).split("=")[1]
            geneName = (row['attributes'].split(";")[1]).split("=")[1]
            geneSource = "Aceview"
            geneBiotype = "protein_coding" ## Aceview assumes peptide sequences can be derived from all transcripts
            geneAttribute="gene_id \"" + geneID + "\"; gene_name \"" + geneName + "\"; gene_source \"Aceview\"; gene_biotype\"protein_coding\""
            attributeList.append(geneAttribute)
            parentList.append(geneID)
            idList.append(geneID)
        if row['feature_type'] == "mRNA":  # If entry is a transcript
            chrList.append(row['chr'])
            sourceList.append(row['source'])
            typeList.append("transcript")
            startList.append(row['start'])
            stopList.append(row['stop'])
            scoreList.append(row['score'])
            strandList.append(row['strand'])
            frameList.append(row['frame'])
            # Reformat attributes -- from GFF3, there is transcript ID, name and parent gene only
            xsID = (row['attributes'].split(";")[0]).split("=")[1]
            xsName = (row['attributes'].split(";")[1]).split("=")[1]
            geneID = (row['attributes'].split(";")[2]).split("=")[1]
            geneName = (row['attributes'].split(";")[2]).split("=")[1]
            xsSource = "Aceview"
            xsBiotype = "protein_coding"  ## Aceview assumes peptide sequences can be derived from all transcripts
            xsAttribute = "gene_id \"" + str(geneID) + "\"; transcript_id \"" + str(xsID) + "\"; gene_name \""\
                          + str(geneName) + "\"; transcript_name \"" + xsName + "\"; transcript_source \""\
                          + xsSource + "\"; gene_biotype \"" + xsBiotype + "\"; transcript_biotype \"" + xsBiotype + "\""
            attributeList.append(xsAttribute)
            parentList.append(geneID)
            idList.append(xsID)
        if row['feature_type'] == "CDS":  # If entry is CDS
            chrList.append(row['chr'])
            sourceList.append(row['source'])
            typeList.append(row['feature_type'])
            startList.append(row['start'])
            stopList.append(row['stop'])
            scoreList.append(row['score'])
            strandList.append(row['strand'])
            frameList.append(row['frame'])
            # Reformat attributes -- from GFF3, there is transcript ID, name and parent gene only
            cdsID = (row['attributes'].split(";")[0]).split("=")[1]
            xsID = (row['attributes'].split(";")[1]).split("=")[1]
            xsName = (row['attributes'].split(";")[1]).split("=")[1]
            proteinID = (row['attributes'].split(";")[1]).split("=")[1]
            geneID = ((row['attributes'].split(";")[1]).split("=")[1]).split('.')[0]
            geneName = ((row['attributes'].split(";")[1]).split("=")[1]).split('.')[0]
            geneSource = "Aceview"
            xsSource = "Aceview"
            xsBiotype = "protein_coding"
            xsAttribute = "ccds_id \"" + cdsID + "\"; gene_id \"" + geneID + "\"; gene_name \"" + geneName +\
                          "\"; transcript_id \"" + xsID + "\"; transcript_name \"" + xsName + "\"; gene_source \"" + \
                          geneSource + "\"; transcript_source \"" + xsSource + "\"; transcript_biotype \"" + \
                          xsBiotype + "\"; protein_id \"" + proteinID + "\";"
            attributeList.append(xsAttribute)
            parentList.append(xsID)
            idList.append(cdsID)
        if row['feature_type'] == "exon":  # If entry is an exon
            exonChr=row['chr']
            exonSource=row['source']
            exonType=row['feature_type']
            exonStart=row['start']
            exonStop=row['stop']
            exonScore=row['score']
            exonStrand=row['strand']
            exonFrame=row['frame']
            exonID = (row['attributes'].split(";")[0]).split("=")[1]
            xsSource = "Aceview"
            geneSource = "Aceview"
            geneBiotype = "protein_coding"
            xsList = ((row['attributes'].split(";")[1]).split("=")[1]).split(",")
            for xs in range(0,len(xsList)):
                geneID=xsList[xs]
                geneName=xsList[xs]
                xsName=xsList[xs]
                xsID=xsList[xs]
                xsAttribute = "exon_id \"" + exonID + "\"; transcript_id \"" + xsID + "\"; transcript_name \"" \
                              + xsName + "\"; gene_id \"" + geneID + "\"; gene_name \"" + geneName\
                              + "\"; transcript_source \"" + xsSource + "\"; gene_source \"" + geneSource \
                              + "\"; gene_biotype \"" + geneBiotype + "\"; transcript_biotype \""\
                              + geneBiotype + "\""
                chrList.append(exonChr)
                sourceList.append(exonSource)
                typeList.append(exonType)
                startList.append(exonStart)
                stopList.append(exonStop)
                scoreList.append(exonScore)
                strandList.append(exonStrand)
                frameList.append(exonFrame)
                attributeList.append(xsAttribute)
                parentList.append(xsID)
                idList.append(exonID)

    gtfDF = pd.DataFrame({'parent_id':parentList, 'feature_id':idList, 'chr': chrList,'source': sourceList, 'feature_type': typeList,
                          'start': startList, 'stop': stopList, 'score': scoreList,
                          'strand': strandList, 'frame': frameList, 'attributes':attributeList })
    gtfDF = gtfDF[['chr','source','feature_type','start','stop','score','strand','frame','attributes']]
    print(gtfDF.head(n=50))
    gtfDF.to_csv(outputFile, encoding='utf-8', index=False, sep='\t', header=None, quoting=csv.QUOTE_NONE)