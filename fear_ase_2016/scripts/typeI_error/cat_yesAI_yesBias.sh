#!/bin/bash
#PBS -M jfear@ufl.edu
#PBS -m n
#PBS -r n
#PBS -q bio
#PBS -o /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/bayesian/
#PBS -e /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/bayesian/
#PBS -l walltime=11:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=500mb
#PBS -t 1-200

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
PYPROG=/scratch/lfs/mcintyre/python.git/catTable.py

OUTDIR=$PROJ/typeI_error/output
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

INPUT=$OUTDIR/split/simulated_YESAI_YESBIASRequal1d5Poisson

python $PYPROG -f ${INPUT}_*.csv \
               --odir $OUTDIR \
               --header
