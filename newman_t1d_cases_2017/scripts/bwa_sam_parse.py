#!/usr/bin/env python
import argparse

## This script parses a sam file from BWA-MEM and outputs a log of alignment counts and percentages.

# Parse command line arguments
parser = argparse.ArgumentParser(description='Parse sam file to get alignment counts.')
parser.add_argument('-sam','--sam_file',dest='sam', action='store', required=True, help='A Sam file to parse [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file for alignment log [Required]')
args = parser.parse_args()

flags=list()

# Open sam file and create a list that contains only the second column from the sam file, (the bitwise flags).
with open(args.sam,'r') as sam:
    for line in sam.readlines():
        cols=line.split('\t')
        flags.append(cols[1])


# Count the flags. These flags are based on BWA sam output, may not be the same for other aligners.
# The flags are different for paired data. There is another python script 'bwa_sam_parse_se.py' for single-end alignments.

unaln=flags.count('77') + flags.count('141') + flags.count('181') + flags.count('121') + flags.count('133') + flags.count('117') + flags.count('69') 

aln=flags.count('99') + flags.count('73') + flags.count('185')  + flags.count('147') + flags.count('83') + flags.count('163') + flags.count('97') + flags.count('137') + flags.count('145') + flags.count('81') + flags.count('161')+ flags.count('177') + flags.count('113') + flags.count('65') + flags.count('129')

ambig=flags.count('337') + flags.count('417') + flags.count('369') + flags.count('433') + flags.count('353') + flags.count('401') + flags.count('371')+ flags.count('355') + flags.count('403') + flags.count('419') + flags.count('339') + flags.count('387') + flags.count('385') + flags.count('323') + flags.count('435') + flags.count('321')

total = unaln + aln

# Get percentages

percent_aln = float (aln) / (total) * 100
percent_unaln = float (unaln) / (total) * 100
percent_ambig = float (ambig) / (total) * 100


# Write the counts to the output.

with open(args.out,'w') as dataout:
    dataout.write('Total reads '+str(total)+'\nAligned '+str(aln)+'\nUnaligned '+str(unaln)+'\nAmbiguous '+str(ambig)+'\nPercent aligned '+str(percent_aln)+'\nPercent unaligned '+str(percent_unaln)+'\nPercent ambiguous '+str(percent_ambig))

