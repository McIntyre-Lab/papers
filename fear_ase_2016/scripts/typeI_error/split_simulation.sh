#!/bin/bash - 
#===============================================================================
#   DESCRIPTION: Take simulations and split them for running on the HPC.
# 
#===============================================================================

set -o nounset                              # Treat unset variables as an error

# Set paths depending on if I am running locally or on HPC
if [[ $HOSTNAME == *ufhpc* ]]; then
    # If on UF HPC
    module load python/2.7.6

    # Set Directories
    PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
else
    # Set Directories
    PROJ=$MCLAB/cegs_ase_paper
fi

# Use split table from mcscript
PYPROG=$PROJ/scripts/mcscript/splitTable.py

for line in r361 r332 r365; do
    OUTDIR=$PROJ/pipeline_output/typeI_error/input/${line}_split
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

    # No AI No Bias
    INPUT=$PROJ/pipeline_output/typeI_error/output/${line}_noai_nobias_sim.csv
    if [ -f $INPUT ]; then
        python $PYPROG -f $INPUT \
                       -o $OUTDIR \
                       --header \
                       --nfiles 500
    fi

    # AI No Bias
    INPUT=$PROJ/pipeline_output/typeI_error/output/${line}_ai_nobias_sim.csv
    if [ -f $INPUT ]; then
    python $PYPROG -f $INPUT \
                   -o $OUTDIR \
                   --header \
                   --nfiles 500
    fi

    # AI Bias
    INPUT=$PROJ/pipeline_output/typeI_error/output/${line}_ai_bias_sim.csv
    if [ -f $INPUT ]; then
    python $PYPROG -f $INPUT \
                   -o $OUTDIR \
                   --header \
                   --nfiles 500
    fi
done
