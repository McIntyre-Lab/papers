#!/bin/bash

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/catTable.py

OUTDIR=$PROJ/typeI_error/output
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

INPUT=$OUTDIR/r101_split/r101_simulated_NOAI_NOBIASRequal1Poisson_results

python $PYPROG -f ${INPUT}_*.csv \
               --odir $OUTDIR \
               --header
