#!/bin/bash - 
#===============================================================================
#   DESCRIPTION: Take misspecification simulations and split them for running
#   on the HPC.
# 
#===============================================================================

set -o nounset                              # Treat unset variables as an error
set -o errexit

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

for line in r361 r332 r365 r101; do
    OUTDIR=$PROJ/pipeline_output/typeI_error/input/${line}_miss_split
    if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

    for QNAME in 35 375 4 425 45 46 47 48 49 50 51 52 53 54 55 575 6 625 65; do
        INPUT=$PROJ/pipeline_output/typeI_error/output/${line}_misspecification_qTrue0${QNAME}_sim.csv
        python $PYPROG -f $INPUT \
                       -o $OUTDIR \
                       --header \
                       --nfiles 500
    done
done
