#!/usr/bin/env python
# Built-in Modules
import os.path
import re
import logging
import argparse
import copy
import sys
import gzip

# Other modules
import vcf as pyvcf

# McIntyre Modules
import mclib
import mclib.vcf2 as mcvcf

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="This script takes a multi-sample vcf and splits the samples into separate vcf files.")
    parser.add_argument("--vcf", dest="vname", action='store', required=True, help="Name of the VCF file, preferably zipped in using bgzip. [Required]")
    parser.add_argument("--outdir", dest="odir", action='store', required=True, help="Name of the output directory. [Required]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log. [Optional]")
    args = parser.parse_args()
    #args = parser.parse_args(['--vcf', '/home/jfear/storage/useful_dmel_data/CEGS.68.lines.raw.SNPs.filtered.set.1.recode.vcf.gz', '--outdir', '/home/jfear/tmp'])
    return(args)

def getName(sample, odir):
    """ Simple function to return a new filename with the sample ID placed in the name """
    if 'w1118' in sample:
        sname = re.sub('w1118_w118', 'w1118',sample)
    elif 'Raleigh' in sample:
        sname = re.sub('Raleigh_','r',sample)
    elif 'Winters' in sample:
        sname = re.sub('Winters_','w',sample)
    else:
        logging.warn('Sample naming does not match cegs please change')
    return os.path.join(odir, sname + '.vcf')

if __name__ == '__main__':
    args = getOptions()
    mclib.logger.set_logger(args.log)

    # Get list of samples
    vcf = mcvcf.Vcf(args.vname)
    samples = vcf.vcf_reader.samples
    numSamples = len(samples)

    for sample in samples:
        sname = getName(sample, args.odir)
        sindex = 0
        with open(sname, 'w') as OUT:
            with gzip.open(args.vname, 'rb') as VCF:
                for row in VCF:
                    if row.startswith("##"):
                        OUT.write(row)
                    elif row.startswith("#"):
                        cols = row.strip().split('\t')
                        sindex = cols.index(sample)
                        myOut = [str(x) for x in cols[0:9] + [cols[sindex]]]
                        OUT.write('\t'.join(myOut) + "\n")
                    else:
                        cols = row.strip().split('\t')
                        myOut = [str(x) for x in cols[0:9] + [cols[sindex]]]
                        OUT.write('\t'.join(myOut) + "\n")
