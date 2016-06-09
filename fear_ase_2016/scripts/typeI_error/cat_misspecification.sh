#!/bin/bash

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/catTable.py

OUTDIR=$PROJ/pipeline_output/typeI_error
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

for LINE in r332 r361 r365 r101
do
    INPUT=$OUTDIR/output/${LINE}_miss_split/${LINE}_misspecification_results_qTrue
    python $PYPROG -f ${INPUT}*.csv \
                   --odir $OUTDIR \
                   --oname ${LINE}_misspecification_results_sim_summary.csv \
                   --header

    sed -i 's/NA//g' $OUTDIR/${LINE}_misspecification_results_sim_summary.csv
done
