#!/usr/bin/env python
import os
import glob

import numpy as np
import scipy.stats as sp

import mclib_Python.gff as mcgff
import mclib_Python.bam as mcbam

mclab = os.getenv('MCLAB')

################################################################################
# Get Annotation information for CR43653
################################################################################
gname = os.path.join('/home/jfear/storage/useful_dmel_data/dmel-all-no-analysis-r5.51.gff')
fly = mcgff.FlyGff(gname)


################################################################################
# Get list of lines and prep output
################################################################################
# Pull a list of lines from the file names
design_file = os.path.join(mclab, 'cegs_sem_sd_paper/design_file/cegsV_line_list.csv')
with open(design_file, 'r') as FH:
    lines = [i.rstrip() for i in FH.readlines()]

# Open output file
oname = os.path.join(mclab, 'cegs_sem_sd_paper/analysis_output/splicing/InR_ncRNA_counts.csv')
OUT = open(oname, 'w')

# Write header
OUT.write(','.join(['sample', 'exon_id', 'mapped_reads', 'read_length', 'region_length', 'region_depth', 'reads_in_region', 'apn', 'rpkm']) + '\n')


################################################################################
# For each line grab pile up information
################################################################################
# Iterate through lines and calculate coverage for the ncRNA exons
for line in lines:
    files = glob.glob('/mnt/storage/cegs_aln/combined/{0}*V*.bam'.format(line))
    summary = list()
    for cfile in files:
        bam = mcbam.Bam(cfile)
        fname = os.path.basename(cfile).split('.')[0]

        # Total number of mapped reads
        numMappedReads = 0
        for read in bam.bam.fetch():
            numMappedReads += 1

        for gene in ('CR43653', 'CR44034'):
            # Look up fbgn
            fbgn = fly.symbol2fbgn[gene]

            # Grab current gene and its atributes
            for exon in fly.db.children(fbgn, featuretype='exon'):
                ename = exon.id
                chrom = exon.chrom
                start = exon.start
                end = exon.end
                length = end - start + 1

                # Pileup in gene region of interest
                pileups = bam.get_pileup(chrom, start, end)
                pileArray = np.array(pileups.values())
                numReadsInRegion = len(list(bam.bam.fetch(chrom, start, end)))

                if len(pileArray) > 0:
                    # Calculate various statistics
                    depth = np.sum(pileArray)
                    apn = depth / float(length)
                    rpkm = (10**9 * numReadsInRegion) / float(numMappedReads * length)
                    #std = np.std(pileArray)
                    #cv = sp.stats.variation(pileArray)
                else:
                    depth = 0
                    apn = 0
                    rpkm = 0
                    std = 0
                    cv = 0

                OUT.write(','.join([str(i) for i in [fname, ename, numMappedReads, '95', length, depth, numReadsInRegion, apn, rpkm]]) + '\n')
