#!/usr/bin/env python

###############################################################################
#
#  NAME: subset_fastq.py
#  DESCRIPTION: This script extracts requested reads from a FASTQ file and
#  outputs them to a new FASTQ file. The default input is a SAM or BAM file,
#  however any text file with the read ID in the first column is acceptable.
#   
#  AUTHOR: Jeremy R. B. Newman
#
################################################################################

# Load packages
import os
import argparse
from Bio import SeqIO

def getOptions():
        """Function to pull in arguments"""
        parser = argparse.ArgumentParser(description='Parse input file to extract desired FASTQ reads.')
        parser.add_argument('--input', dest='input', action='store', required=True, help='Input file, usually SAM/BAM [Required]')
        parser.add_argument('--in-fq', dest='inFQ', action='store', required=True, help='Input FASTQ file [Required]')
        parser.add_argument('--out-fq', dest='outFQ', action='store', required=True, help='Output FASTQ file [Required]')
        args = parser.parse_args()
        return (args)

# Get arguments
args = getOptions()

# Import data
samFile=[read.strip().split() for read in open(args.input).readlines()]

# Get list of reads to extract
wantReads=[]
for seq in range(0,len(samFile)):
    wantReads.append(samFile[seq][0])

# Open original FASTQ
records = list(SeqIO.parse(args.inFQ, "fastq"))

# Extract reads
testout=[]
for read in range(0,len(wantReads)):
    for fq in range(0,len(records)):
        if ''.join(wantReads[read]) == records[fq].id:
            testout.append(records[fq])

# Write output
SeqIO.write(testout, args.outFQ, "fastq")

