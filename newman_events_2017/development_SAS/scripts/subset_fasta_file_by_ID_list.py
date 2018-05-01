#!/usr/bin/env python

#
#    DESCRIPTION: This script takes a fASTA file and subsets it using a list of Fasta IDs.
#    It creates a new FASTA file containing only the reads in the subset list.
#
#    AUTHOR: Alison Morse
#

import os,csv,sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Parse command line arguments
parser = argparse.ArgumentParser(description='Subset FA file using a list of FA IDs')
parser.add_argument("-i", "--input_file", dest="input_file", required=True, help="fa file to subset from")
parser.add_argument("-l", "--id_file", dest="id_file", required=True, help="list if IDs want to subset")
parser.add_argument("-o", "--output_file", dest="output_file", required=True, help="output file for the subsetted fa reads")
args = parser.parse_args()

wanted = set()
with open(args.id_file, "r") as id:
    for line in id:
        line = line.strip()
        if line !="":
            wanted.add(line)

fa_sequences = SeqIO.parse(open(args.input_file, "r"), "fasta")
with open(args.output_file, "w") as out:
    for seq in fa_sequences:
        if seq.id in wanted:
            SeqIO.write([seq], out, "fasta")
