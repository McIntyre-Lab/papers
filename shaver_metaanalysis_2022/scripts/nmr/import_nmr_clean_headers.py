#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:29:17 2020

@author: ammorse
"""

import argparse
import os
import pandas as pd
from pandas import DataFrame

## import and rename headers in Carters feature list and compound list so can match 

pathInOut="/home/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/data_from_group/nmr_sweet16"

## import dataset 
df = pd.read_csv(os.path.join(pathInOut, "ph_NMR_CDCL3_14Aug2020_flagsdropped.txt"), sep='\t', header=0)
list(df)

# replace CDCL3 with nothing
df2 = df.rename(columns=lambda x: x.replace('NMR_CDCL3_', ''))


outS = 'ph_NMR_CDCL3_shortNames.tsv'
outFileS = os.path.join(pathInOut, outS)
df2.to_csv(outFileS, index=False, sep='\t')


## import dataset 
df3 = pd.read_csv(os.path.join(pathInOut, "ph_NMR_D2O_14Aug2020.txt"), sep='\t', header=0)
list(df3)

# replace CDCL3 with nothing
df4 = df3.rename(columns=lambda x: x.replace('NMR_D2O_', ''))


outS2 = 'ph_NMR_D2O_shortNames.tsv'
outFileS2 = os.path.join(pathInOut, outS2)
df4.to_csv(outFileS2, index=False, sep='\t')