#!/bin/bash
#PBS -M jrbnewman@ufl.edu
#PBS -m n
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/sugrue/scripts/PBS_LOGS/sas
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=16G
#PBS -t 1

## 25
module load sas/9.3

#Set up variables

PROJ=/scratch/lfs/sugrue
SCRIPTS=$PROJ/scripts/sas
LOGS=$PROJ/sas_analysis/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_analysis/sas_data
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi

#PRG0=import_and_format_junctions.sas
#PRG0=junctions_flag_counts.sas
#PRG0=junctions_transpose_data.sas
#PRG0=junctions_calc_means_jrbn2.sas
#PRG0=junction_assemble_data.sas
#PRG0=junctions_anova.sas
#PRG0=junction_make_final_results_jrbn2.sas
PRG0=junctions_fdr_final.sas





TMPDIR=/scratch/lfs/sugrue/sas_temp/splicing_${PBS_ARRAYID}

if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi

cd $SASDATA

sas -log $LOGS/${PRG0}_${PBS_ARRAYID}.log -work $TMPDIR -sysin $SCRIPTS/$PRG0 -print $LOGS/${PRG0}_${PBS_ARRAYID}.lst -sysparm ${PBS_ARRAYID}

rm $TMPDIR/*


