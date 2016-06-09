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

INPUT=$PROJ/emp_bayesian/PG_model/split
OUTDIR=$PROJ/emp_bayesian/PG_model
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi 

# combine table
python $PYGIT/catTable.py -f $INPUT/split_*.csv --odir $OUTDIR --oname PG_emp_bayesian_results.csv --header 
