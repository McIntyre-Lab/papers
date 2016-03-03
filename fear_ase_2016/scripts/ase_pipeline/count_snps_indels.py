#!/usr/bin/env python
# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys
import re

# Add-on packages

# McLab Packages
import mclib


def getOptions():
    """ Function to pull in arguments """

    description = """ This script takes a VCF 4.2 file pulls out snps/indels and depth into a csv """
    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    group1 = parser.add_argument_group(description="Input Files")
    group1.add_argument("--vcf", dest="vname", action='store', required=False, help="Name of the uncompressed VCF file v4.2. [stdin]")

    group2 = parser.add_argument_group(description="Output Files")
    group2.add_argument("-o", dest="oname", action='store', required=False, help="Name of output csv. [stdout]")
    group2.add_argument("--log", dest="log", action='store', required=False, help="Name of the LOG file [stderr]") 

    parser.add_argument("--debug", dest="debug", action='store_true', required=False, help="Enable debug output.") 

    args = parser.parse_args()
    return(args)

def inputOutput(args):
    """ If an input/output files were given then use them, else use
    stdin/stdout 
    
    Arguments:
    ----------
    args (obj): command line arguments.

    Returns:
    --------
    FH = a input file handlder
    OUT = a output file handler
    """

    if args.vname:
        logger.info('Reading file %s' % args.vname)
        FH = open(args.vname, 'r')
    else:
        logger.info('Reading STDIN')
        FH = sys.stdin
    if args.oname:
        logger.info('Writing to %s' % args.vname)
        OUT = open(args.oname, 'w')
    else:
        logger.info('Writing to STDOUT')
        OUT = sys.stdout
    return FH, OUT

def parseLine(line):
    """ Take a VCFv2.4 row and split in up into pieces """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.rstrip().split('\t')

    # split info into dictionary
    infoDict = dict()
    for item in info.split(';'):
        if item == 'INDEL':
            # If it is an indel then set key and value
            key = 'INDEL'
            value = '1';
        else:
            # Split up based on '='
            key, value = item.split('=')

        # If value is a list then split it up. Figure out if
        # these are numbers or strings
        if ',' in value:
            try:
                val = [int(x) for x in value.split(',')]
            except:
                val = value.split(',')
        else:
            try: 
                val = int(value)
            except:
                val = value
        infoDict[key] = val

    # Sometimes ref and alt have multiple alleles, split them into a list
    refList = ref.split(',')
    altList = alt.split(',')

    # Samtools includes the value <X> for all SNPs and it includes a count in
    # the DPR field. I want to remove these.
    if '<X>' in altList:
        altList.remove('<X>')
        infoDict['DPR'].pop()

    # Merge ref and alt back together, this time using a | as separator so I
    # don't interfer with the csv format
    refClean = '|'.join(refList)
    altClean = '|'.join(altList)

    return chrom, pos, refClean, altClean, infoDict

def main(args):
    # Get Input/Output file handlders
    FH, OUT = inputOutput(args)

    header = ','.join(['chrom', 'pos', 'ref', 'alt', 'totalCount', 'refCount', 'altCount', 'flagIndel']) + '\n'
    OUT.write(header)

    # Make sure the input vcf file is v4.2
    logger.info('Checking VCF version')
    line1 = FH.next()
    vcfVersion = line1.split('=')[1].rstrip()
    if vcfVersion != 'VCFv4.2':
        logger.error('This script was written only for vcf version 4.2')
        sys.exit(1)
    else:
        logger.info('Your VCF is version 4.2')

    # Iterate over lines and grab usefule info
    logger.info('Parsing VCF file')
    for line in FH:
        # Skip header lines
        if line.startswith('#'):
            continue

        # split the line into its parts
        chrom, pos, ref, alt, infoDict = parseLine(line)
        refCount = infoDict['DPR'][0]
        altCount = '|'.join([str(x) for x in infoDict['DPR'][1:]])
        totalCount = sum(infoDict['DPR'])

        if 'INDEL' in infoDict:
            flagIndel = 1
        else:
            flagIndel = 0

        # Write output as csv
        if alt:
            myOut = ','.join([str(x) for x in [chrom, pos, ref, alt, totalCount, refCount, altCount, flagIndel]]) + '\n'
            OUT.write(myOut)

    FH.close()
    OUT.close()


if __name__ == '__main__':
    # Turn on Logging if option -g was given
    args = getOptions()

    # Turn on logging
    logger = logging.getLogger()
    if args.debug:
        mclib.logger.setLogger(logger, args.log, 'debug')
    else:
        mclib.logger.setLogger(logger, args.log)

    # Run Main part of the script
    main(args)
    logger.info("Script complete.")
