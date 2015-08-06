#!/usr/bin/env python

import numpy as np
import pandas as pd


def build_allele_dict():
    """ Take a sheet and build a dictionary with:
    [gene][allele] = count
    """




fname = '/home/jfear/mclab/cegs_sem_sd_paper/from_matt/DSRP_and_CEGS_haps_1-6-15.xlsx'
data = pd.ExcelFile(fname)

dspr = data.parse('DSRP_haps')
f1 = data.parse('CEGS_haps')
data
