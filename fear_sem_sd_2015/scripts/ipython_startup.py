from __future__ import division
# Run these command when starting an ipython notebook. This is simply for
# convience, so you don't have to do all of this in every python notebook.

# Standard Import
print("Importing commonly used libraries: os, sys, numpy as np, scipy as sp, pandas as pd, matplotlib as mp, matplotlib.pyplot as plt, datetime as dt")
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
import datetime as dt

# Set current date
TODAY = dt.date.today().strftime("%Y%m%d")
# Set up paths
## Get mclab env
### I have an environmental variable MCLAB that points to a share directory that
### I use across all of my machines.
MCLAB = os.getenv('MCLAB')

## Get Current project directory
### This is the project folder on the shared drive.
PROJ = os.path.join(MCLAB, 'cegs_sem_sd_paper')

print("Creating project level variables: MCLAB = {}, PROJ = {}, TODAY = {}".format(MCLAB, PROJ, TODAY))

## Add the mclib_Python libraries to PYTHONPATH
### You can point to any copy of this library. I like to keep a copy for each
### project so I know exactly what has been run.
sys.path.append(os.path.join(PROJ, 'scripts/mclib_Python'))

## Add project level library

# Run IPython Magics
from IPython import get_ipython
ipython = get_ipython()

## Have matplotlib plot inline
ipython.magic("matplotlib inline")

## Turn on autoreload
### This allows you to modify functions and have them automatically reloaded
### when they are called.
ipython.magic("load_ext autoreload")
ipython.magic("autoreload 2")
