GENENAME = 'ps'
import matplotlib
matplotlib.use('Agg')

import mclib_Python
from mclib_Python import gff as mcgff
from mclib_Python import bam as mcbam
import numpy as np
from glob import glob
import pandas as pd
import seaborn
import matplotlib.pyplot as plt

# Import GFF File and Get gene location
db = mcgff.FlyGff('/scratch/lfs/mcintyre/references/dmel_fb551/dmel-all-no-analysis-r5.51.gff')
gene = mcgff.FlyGene(GENENAME, db)

def get_pileup(fname, chrom, start, end):
    """ Function to pull out reads from BAM files. """
    bam = mcbam.Bam(fname)
    pileup = bam.get_pileup(chrom, start, end)
    return pd.Series(pileup)

# Get list of lines
with open('/scratch/lfs/mcintyre/cegs_ase_paper/design_files/CEGS_list_68_lines.txt', 'r') as FH:
    lines = FH.read().rstrip('\n').split('\n')

# Iterate over lines
matrix = list()
for LINE in lines:    
    pileups = list()
    for FILE in glob('/scratch/lfs/mcintyre/cegs_oe/bam_fb551_genome_nodup/{}_*.sorted.bam'.format(LINE)):
        pileups.append(get_pileup(FILE, gene.chrom, gene.start, gene.end))

    # Sum coverage
    matrix.append(pd.concat(pileups, axis=1).fillna(0).sum(axis=1))

df = pd.concat(matrix, axis=1).fillna(0)

df.columns = lines

fig, ax = plt.subplots(1, 1, figsize=(20, 10))
seaborn.heatmap(df.T, ax=ax)
fig.savefig('/scratch/lfs/mcintyre/cegs_ase_paper/wiggles/ps_heatmap.png')
