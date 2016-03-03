#!/bin/bash

#PBS -M jrbnewman@ufl.edu
#PBS -m n
#PBS -w newgroup=concannon
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/patcon/fnew/scripts/PBS_LOGS/sas
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=8G



####PBS -t 1-5,7-54

#### 50
module load sas/9.3

#Set up variables

PROJ=/scratch/lfs/patcon/jnewman
SCRIPTS=$PROJ/scripts/sas
LOGS=$PROJ/sas_analysis/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_analysis/sas_data
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi

#PRG0=export_all_files.sas
#PRG0=merge_${PBS_ARRAYID}.sas
#PRG0=add_covariates.sas
#PRG0=pull_event_info.sas
PRG0=analysis_sum_tech_reps.sas
#PRG0=analysis_design_file.sas
#PRG0=analysis_import_coverage_counts.sas
#PRG0=analysis_import_coverage_counts_splicing.sas
#PRG0=analysis_sum_tech_reps_splicing.sas
#PRG1=analysis_quartile_normalization.sas
#PRG2=analysis_flag_isoform_on_off.sas
#PRG3=analysis_exp_means.sas
#PRG4=analysis_model_all_on_apn0.sas
#PRG5=analysis_flag_fdrs.sas
#PRG6=analysis_create_final_dataset.sas

TMPDIR=/scratch/lfs/patcon/jnewman/sas_temp
if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi

cd $SASDATA

#sas -log $LOGS/$PRG0_${PBS_ARRAYID}.log -work $TMPDIR -sysin $SCRIPTS/$PRG0 -sysparm ${PBS_ARRAYID}
sas -log $LOGS/$PRG0.log -work $TMPDIR -sysin $SCRIPTS/$PRG0
#sas -log $LOGS/$PRG1.log -work $TMPDIR -sysin $SCRIPTS/$PRG1
#sas -log $LOGS/$PRG2.log -work $TMPDIR -sysin $SCRIPTS/$PRG2
#sas -log $LOGS/$PRG3.log -work $TMPDIR -sysin $SCRIPTS/$PRG3
#sas -log $LOGS/$PRG4.log -work $TMPDIR -sysin $SCRIPTS/$PRG4
#sas -log $LOGS/$PRG5.log -work $TMPDIR -sysin $SCRIPTS/$PRG5
#sas -log $LOGS/$PRG6.log -work $TMPDIR -sysin $SCRIPTS/$PRG6

rm $TMPDIR/*

