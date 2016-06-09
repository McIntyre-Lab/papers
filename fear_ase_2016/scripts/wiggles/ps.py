GENENAME = 'ps'
import matplotlib
matplotlib.use('Agg')

import mclib
from mclib import gff as mcgff
from mclib import bam as mcbam
import numpy as np
from glob import glob
import pandas as pd
import seaborn
import matplotlib.pyplot as plt

# Import GFF File and Get gene location
db = mcgff.FlyGff('/home/agerken/mclab/useful_dmel_data/flybase551/flybase_files/dmel-all-no-analysis-r5.51.gff')
gene = mcgff.FlyGene(GENENAME, db)

def get_pileup(fname, chrom, start, end):
    """ Function to pull out reads from BAM files. """
    bam = mcbam.Bam(fname)
    pileup = bam.get_pileup(chrom, start, end)
    return pd.Series(pileup)

# Get list of lines
with open('/home/jfear/lines.txt', 'r') as FH:
    lines = FH.read().rstrip('\n').split('\n')

# Iterate over lines
matrix = list()
for LINE in lines:    
    pileups = list()
    for FILE in glob('/media/agerken/backup/cegs_oe/combined/{}_*.sorted.bam'.format(LINE)):
        pileups.append(get_pileup(FILE, gene.chrom, gene.start, gene.end))

    # Sum coverage
    matrix.append(pd.concat(pileups, axis=1).fillna(0).sum(axis=1))

df = pd.concat(matrix, axis=1).fillna(0)

df.columns = lines

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
seaborn.heatmap(df.T, ax=ax)
fig.savefig('/home/agerken/cegs_ase_paper/pipeline_output/wiggle_plots/ps_heatmap.png')
