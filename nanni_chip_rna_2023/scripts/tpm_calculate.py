#!/usr/bin/env python

import argparse
import logging
import csv
import collections
import numpy as np
from argparse import RawDescriptionHelpFormatter
csv.field_size_limit(1000000000)

def getOptions():
    """ Function to pull in arguments """

    description="""This script can be used to calculates coverage (APN and TPM) in two different settings:\n
    (1) Coverage can be calculated across an entire genomic region. To do
        this a 3-column bed file must be provided (Try fasta2bed.py). 
        col[0] = chromosome/featureID name (eg., chr2L or S7_SI)
        col[1] = start position (i.e., '0') 
        col[2] = end position (i.e., length) 

    (2) Coverage can also be calculated by excising specific exons/features
        from a genome. For example if you have aligned to the genome, but want
        coverage at the exon level. For this a 4-column bed must be provided.
        col[0] = chromosome name (eg., chr2L)
        col[1] = exon/feature start position (eg., 2929)
        col[2] = exon/feature end position (eg., 3090)
        col[3] = exon/feature name (eg., S7_SI) 
        IMPORTANT: Setting 2 requires a lot of RAM ~10-12GB for calculating
        coverage using fly feature """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--mpileup", dest="mname", action='store', required=True, help="mpileup file [Required]",metavar="MPILEUP_FILE")
    parser.add_argument("-s", "--sam", dest="sname", action='store', required=True, help="sam alignment file [Required]", metavar="SAM_FILE")
    parser.add_argument("-b", "--bed", dest="bname", action='store', required=True, help="bed file (3 or 4 columns) [Required]", metavar="BED_FILE")
    parser.add_argument("-c", "--cv", dest="cv", action='store_true', required=False, help="A flag to indicate if you want output for the coefficient of variation [OPTIONAL]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    parser.add_argument("-o", "--out", dest="out", action='store', required=True, help="Out File", metavar="OUT_FILE")
    args = parser.parse_args()
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

# SAM Functions
def read_sam(args):
    """ Read SAM file to get read length and number of mapped reads. Note: if
        you have allowed ambiguous mapping then reads are counted multiple times.  """
    logging.info("Reading the SAM File '%s'." % args.sname)
    num_mapped_reads = 0
    read_length = 0
    with open(args.sname,'r') as SAM:
        for row in SAM:
            if not row.startswith('@'):
                record = row.strip().split('\t')
                if record[1] != 4:          # only look at aligned reads, a 4 in column 2 of a single end aligned sam file indicates an unaligned read (see sam format)
                    num_mapped_reads += 1
                read_length = max(read_length,len(record[9]))    # find the maximum read length
    return(num_mapped_reads,read_length)

# BED Functions
def read_bed(args):
    """ Read BED file and create a dictionary containing all information """
    logging.info("Reading the BED File '%s'." % args.bname)
    bdict = collections.defaultdict(dict)
    with open(args.bname,'r') as BED:
        reader = csv.reader(BED,delimiter='\t')
        for row in reader:
            if len(row) == 4:                    # If BED file has 4 columns
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                length = end - start
                feature = row[3]
            elif len(row) == 3:                  # If BED file has 3 columns
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                length = end
                feature = row[0]
            else:
                logging.error("I can only handle 3 or 4 column bed files. See Help for descriptions")
                exit(-1)
            bdict[feature]['chrom'] = chrom
            bdict[feature]['start'] = start
            bdict[feature]['end'] = end
            bdict[feature]['region_length'] = length + 1  # convert back to 1 based scale
            bdict[feature]['count'] = np.zeros(length)    # create a holding vector of 0's as long as the region, I will replace the 0's with counts from the mpileup
    cdict = collections.defaultdict(dict)
    for feature in bdict:
        chrom = bdict[feature]['chrom']
        start = bdict[feature]['start']
        end = bdict[feature]['end']
        cdict[chrom].update(dict((x,feature) for x in xrange(start,end+1))) # create a look up dictionary by chromosome. This will make parsing the mpileup faster.
    return(bdict,cdict)

# MPILEUP Functions
def read_mpileup(args,bdict,cdict):
    """ Read mpileup and store depth and length into dictionary """
    logging.info("Reading the Mpileup File '%s'." % args.mname)
    with open(args.mname, 'r') as MPILEUP:
        reader = csv.reader(MPILEUP, delimiter='\t',quoting=csv.QUOTE_NONE)
        for row in reader:
            mchrom = row[0]
            mpos = int(row[1]) - 1  # mpileups are 1-based
            mdepth = int(row[3])
            try:
                feature = cdict[mchrom][mpos]
                loc = mpos - bdict[feature]['start']
                bdict[feature]['count'][loc] = mdepth
            except:
                continue

# Coverage Functions
def calc_coverage(bdict,num_mapped_reads,read_length):
    """ Calculate different coverage metrics: Estimate number of reads in
    region, Average per nucleotide coverage (APN), average coverage across region (mean),
    standard deviation of coverage in region (std), and coefficient of
    variation (cv), reads per kilobase(rpk), per million scaling factor (per_million_scale), and transcripts per kilobase million (tpm). """
    logging.info("Calculating Coverage Counts")
    rpk_sum=0.0
    for feature in bdict:
        depth = np.sum(bdict[feature]['count'])
        bdict[feature]['depth'] = int(depth)
        bdict[feature]['mean'] = np.mean(bdict[feature]['count'])
        bdict[feature]['std'] = np.std(bdict[feature]['count'])
        if depth != 0:
            bdict[feature]['reads_in_region'] = depth / float(read_length)  # Estimate the number of reads in region based on depth/read_length. Multiplying by 1.0 to tell python to use decimals.
            bdict[feature]['apn'] = depth / float(bdict[feature]['region_length'])  # Calculate average per nucleotide coverage APN (depth in region / length of region). Multiplying by 1.0 to tell python to use decimals.
            bdict[feature]['cv'] = bdict[feature]['std'] / bdict[feature]['mean'] # Calculate the coefficient of variation
            bdict[feature]['rpk'] = bdict[feature]['reads_in_region'] / (float(bdict[feature]['region_length']) / 1000) # Calculate reads per kilobase (Divide read count by region length divided by 1000)
            rpk_sum= rpk_sum + bdict[feature]['rpk']
        else:
            # if there is no coverage set everything to 0
            bdict[feature]['reads_in_region'] = 0
            bdict[feature]['apn'] = 0
            bdict[feature]['cv'] = 0
            bdict[feature]['rpk'] = 0
            
    for feature in bdict:
        if depth != 0:
            bdict[feature]['per_million_scale'] = rpk_sum / 1000000 # Calculate per million scaling factor using the sum of all rpk values (rpk/1000000)
            bdict[feature]['tpm'] = bdict[feature]['rpk'] / bdict[feature]['per_million_scale'] # Calculate transcripts per kilobase million (rpk / per million scaling factor)
        else:
            bdict[feature]['per_million_scale'] = 0
            bdict[feature]['tpm'] = 0
            
# Output Functions
def writeOutput(args,bdict,num_mapped_reads,read_length):
    """ I tried writing output using the CSV module, but this did not behave
        well with SAS downstream. So I opted for the brute force method. """
    logging.info("Writing Output")

    if args.cv:
        header = ['featureID','mapped_reads','read_length','region_length','region_depth','reads_in_region','apn','mean','std','cv','tpm']
        with open(args.out, 'wb') as OUT:
            OUT.write(','.join(header) + '\n')
            for key in bdict:
                OUT.write(','.join(str(x) for x in [key,num_mapped_reads,read_length,bdict[key]['region_length'],bdict[key]['depth'],bdict[key]['reads_in_region'],bdict[key]['apn'],bdict[key]['tpm'],bdict[key]['mean'],bdict[key]['std'],bdict[key]['cv']]) + '\n')
    else:
        header = ['featureID','mapped_reads','read_length','region_length','region_depth','reads_in_region','apn','tpm']
        with open(args.out, 'wb') as OUT:
            OUT.write(','.join(header) + '\n')
            for key in bdict:
                OUT.write(','.join(str(x) for x in [key,num_mapped_reads,read_length,bdict[key]['region_length'],bdict[key]['depth'],bdict[key]['reads_in_region'],bdict[key]['apn'],bdict[key]['tpm']]) + '\n')


def main():
    """ MAIN Function to execute everything """
    args = getOptions()
    if args.log:                                         # Turn on Logging if option -g was given
        setLogger(args.log,logging.INFO)

    num_mapped_reads, read_length = read_sam(args)       # Use SAM file to count the number of mapped reads and the max read length
    bdict,cdict = read_bed(args)                               # Read through BED file and create dictionary to sort all information.
    read_mpileup(args,bdict,cdict)                             # Read Mpileup file and populate the bdict with pileup information
    calc_coverage(bdict,num_mapped_reads,read_length)    # Use information stored in bdict to calculate coverage (APN and TPM) and other measures for the genomic region
    writeOutput(args,bdict,num_mapped_reads,read_length) # Write output to CSV file

if __name__=='__main__':
    main()
    logging.info("Script Complete")