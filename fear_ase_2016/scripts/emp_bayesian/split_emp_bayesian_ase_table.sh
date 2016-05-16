#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: Split the ASE table into chunks of rows for faster processing
#   in the Bayesian machine.
# 
#===============================================================================

# Load python
module load python/2.7.6

# Set Variables
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYGIT=/scratch/lfs/mcintyre/python.git

INPUT=$PROJ/emp_bayesian/input/ase_dataset_for_bayesian.csv

OUTDIR=$PROJ/emp_bayesian/input/split
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

# Split table
python $PYGIT/splitTable.py -f $INPUT -o $OUTDIR --prefix split --header --nfiles 500
