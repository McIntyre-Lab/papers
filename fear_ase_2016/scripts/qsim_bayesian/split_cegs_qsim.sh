#!/bin/bash

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/splitTable.py

OUTDIR=$PROJ/qsim/input/split
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

INPUT=$PROJ/qsim/output/ase_dataset_for_bayesian_w_qsim.csv

python $PYPROG -f $INPUT \
               -o $OUTDIR \
               --header \
               --nfiles 500
