#!/bin/bash

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/catTable.py

OUTDIR=$PROJ/qsim/output
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

INPUT=$OUTDIR/split/ase_dataset_for_bayesian_w_qsim

python $PYPROG -f ${INPUT}_*.csv \
               --odir $OUTDIR \
               --header
