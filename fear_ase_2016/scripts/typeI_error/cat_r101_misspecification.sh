#!/bin/bash

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/catTable.py

OUTDIR=$PROJ/typeI_error/output
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

INPUT=$OUTDIR/r101_miss_split/r101_simulateddata_forTIER_FigureNOAIPoisson_q_true_results_0d

python $PYPROG -f ${INPUT}*.csv \
               --odir $OUTDIR \
               --header
