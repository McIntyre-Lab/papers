#!/bin/bash

#PBS -M jrbnewman@ufl.edu
#PBS -m n
#PBS -w newgroup=concannon
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/patcon/fnew/scripts/PBS_LOGS/sas
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=16G
#PBS -t 1


## 11
module load sas/9.3

#Set up variables

PROJ=/scratch/lfs/patcon/jnewman
SCRIPTS=$PROJ/scripts/sas
LOGS=$PROJ/sas_analysis/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_analysis/sas_data2
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi

#PRG0=analysis_import_coverage_counts_splicing2.sas
#PRG0=analysis_splicing_flag_splicing_on_off.sas
#PRG0=analysis_splicing_calc_means.sas
#PRG0=analysis_splicing_run_anova.sas
#PRG0=analysis_splicing_fdr.sas

#PRG0=analysis_splicing_counts_for_eqtls.sas
PRG0=analysis_splicing_export_counts.sas

TMPDIR=/scratch/lfs/patcon/jnewman/sas_temp/splicing_${PBS_ARRAYID}

if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi

cd $SASDATA

sas -log $LOGS/${PRG0}_${PBS_ARRAYID}.log -work $TMPDIR -sysin $SCRIPTS/$PRG0 -sysparm ${PBS_ARRAYID}

rm $TMPDIR/*

