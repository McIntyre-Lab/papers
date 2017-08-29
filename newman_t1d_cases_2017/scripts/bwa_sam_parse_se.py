#!/usr/bin/env python
import argparse

## This script parses a sam file from BWA-MEM and outputs a log of alignment counts and percentages.

# Parse command line arguments
parser = argparse.ArgumentParser(description='Parse sam file to get alignment counts.')
parser.add_argument('-sam','--sam_file',dest='sam', action='store', required=True, help='A Sam file to parse [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for alignment log [Required]')
args = parser.parse_args()

flags=list()

# Open sam file and create a list that contains only the second column from the sam file, the flags.
with open(args.sam,'r') as sam:
    for line in sam.readlines():
        cols=line.split('\t')
        flags.append(cols[1])


# Count the flags. These flags are based on BWA sam output, may not be the same for other aligners.

unaln=flags.count('4')
aln=flags.count('0') + flags.count('16')
ambig=flags.count('256') + flags.count('272')
total = unaln + aln

percent_aln = float (aln) / (total) * 100
percent_unaln = float (unaln) / (total) * 100
percent_ambig = float (ambig) / (total) * 100


# Write the counts to the output log.

with open(args.out,'w') as dataout:
    dataout.write('Total reads '+str(total)+'\nAligned '+str(aln)+'\nUnaligned '+str(unaln)+'\nAmbiguous '+str(ambig)+'\nPercent aligned '+str(percent_aln)+'\nPercent unaligned '+str(percent_unaln)+'\nPercent ambiguous '+str(percent_ambig))

