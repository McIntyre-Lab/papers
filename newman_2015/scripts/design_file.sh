#!/bin/bash

#PBS -M fnew@ufl.edu
#PBS -m n
#PBS -w newgroup=concannon
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/patcon/fnew/scripts/PBS_LOGS/sas
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=10G

module load sas/9.3

#Set up variables

PROJ=/scratch/lfs/patcon/fnew
SCRIPTS=$PROJ/scripts/sas
LOGS=$PROJ/sas_analysis/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_analysis/sas_data
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi


PRG0=analysis_design_file.sas

TMPDIR=/scratch/lfs/fnew/sas_temp
if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi

cd $SASDATA

sas -log $LOGS/$PRG0.log -work $TMPDIR -sysin $SCRIPTS/$PRG0

rm $TMPDIR/*

