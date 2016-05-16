#!/usr/bin/env python
""" Combines all level 2 filtered SNPs/InDels into a single table """

import os
import sys
from glob import glob
import numpy as np
import pandas as pd
import cPickle as pkl

MCLAB = os.getenv('MCLAB')
PROJ = os.path.join(MCLAB, 'cegs_ase_paper')

sys.path.append(os.path.join(PROJ, 'scripts'))

from mclib_Python import vcf2 as mcvcf

# VCF files are stored here
WORK = '/home/jfear/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files'


def main():
    df = pd.DataFrame(columns=['line', 'chrom', 'pos', 'base'])
    for fname in glob(os.path.join(WORK, '*.gz')):
        geno = os.path.basename(fname).split('_')[0]
        if geno != 'w1118':
            vcf = mcvcf.Vcf(fname)
            for row in vcf.vcf_reader:
                chr = row.CHROM
                pos = row.POS
                alt = row.ALT[0]
                ref = row.REF

                if alt == ref:
                    base = np.nan
                else:
                    base = alt

                df.loc[df.shape[0], :] = [geno, chr, pos, base]
    pkl.dump(df, open(os.path.join(PROJ, 'pipeline_output/similarity/variants.pkl', 'wb')))


if __name__ == '__main__':
    main()

