#!/usr/bin/env python

# Take in BBmerge merged fastq output and split reads into >=50% original read length and <50%

from Bio import SeqIO
import argparse
import os
import sys


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Take in BBmerge merged fastq output and split reads into >=50% original read length and <50%")

    # Input arguments
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Merged FASTQ file from BBmerge")
    parser.add_argument("-r", "--readLength", dest="readLength", required=True, type=int, help="Original read length")

    # Output arguments
    parser.add_argument("-g", "--greaterEqualDir", dest="GEdir", required=True, help="Path to directory for the reads that are greater than or equal to 50% the original read length")
    parser.add_argument("-l", "--lessThanDir", dest="Ldir", required=True, help="Path to directory for the reads that are less than 50% the original read length")

    args = parser.parse_args()
    return args

def main():
    # Check that input file and output directories exist
    if not os.path.isfile(args.inFile) or not os.path.isdir(args.GEdir) or not os.path.isdir(args.Ldir):
        raise FileNotFoundError

    # Get filename for output files of the same name in the new output directories
    filename = os.path.basename(args.inFile)

    # Set cutoff length
    cutoff = args.readLength/2
    print("Merged Read cutoff (50% RL) = " + str(cutoff))

    # Import fastq input file and split by sequence length
    listGE = []
    listL = []
    for seqRecord in SeqIO.parse(args.inFile, "fastq"):
        if len(seqRecord) >= cutoff:
            listGE.append(seqRecord)
        else :
            listL.append(seqRecord)
    print("Count of sequences found >= 50% RL = " + str(len(listGE)))
    print("Count of sequences found < 50% RL = " + str(len(listL)))

    SeqIO.write(listGE, args.GEdir + "/" + filename, "fastq")
    SeqIO.write(listL, args.Ldir + "/" + filename, "fastq")

if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()
