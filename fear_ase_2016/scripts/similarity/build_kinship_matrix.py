#!/usr/bin/env python
"""
Kjong has a kinship matrix, but it does not have the labels associated with it.
I take his matrix and associate the labels and output a copy of the matrix for
use in other programs.
"""

# Built-in packages
import os
import cPickle as pkl
import h5py
import numpy as np
import scipy as sp
import pandas as pd

MCLAB = os.getenv('MCLAB')
PROJ = os.path.join(MCLAB, 'cegs_ase_paper')


if __name__ == '__main__':
    kName = os.path.join(PROJ, 'external_data/similarity/kinship.pickle')
    idName = os.path.join(PROJ, 'external_data/similarity/CEGS.216.lines.NO_DPGP4.GATK.SNP.HETS.FILTERED.11-6-13_imputed_v2.hdf5')

    # Load the kinship matrix from the pickle file
    with open(kName, 'rb') as IN:
        kin = pkl.load(IN)

    with h5py.File(idName, 'r') as IN:
        ids = IN['gt_ids'][:]
        ids = [x.replace('Raleigh_', 'r') for x in ids]
        ids = [x.replace('Winters_', 'w') for x in ids]
        ids = [x.replace('w1118_w118', 'w1118') for x in ids]

    df = pd.DataFrame(kin, columns=ids, index=ids)
    df.to_csv(os.path.join(PROJ, 'pipeline_output/similarity/kinship_matrix.csv'))
