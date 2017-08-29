#!/bin/bash

#PBS -M jrbnewman@ufl.edu
#PBS -m n
#PBS -w newgroup=concannon
#PBS -q bio
#PBS -r n
#PBS -j oe
#PBS -o /scratch/lfs/patcon/jnewman/scripts/PBS_LOGS/sas
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=16G
#PBS -t 1

### 25 for eqtls

## 11
module load sas/9.4

#Set up variables

PROJ=/scratch/lfs/patcon/jnewman
SCRIPTS=$PROJ/scripts/sas/eqtls
LOGS=$PROJ/sas_analysis/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_analysis/eqtls
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi

EQTL=$PROJ/eqtls
if [ ! -e $EQTL ]; then mkdir -p $EQTL; fi

RESULTS=$EQTL/results
if [ ! -e $RESULTS ]; then mkdir -p $RESULTS; fi


## Splicing programs
#PRG0=analysis_import_coverage_counts_splicing2.sas
#PRG0=analysis_splicing_flag_splicing_on_off.sas
#PRG0=analysis_splicing_calc_means.sas
#PRG0=analysis_splicing_counts_for_eqtls.sas
#PRG0=analysis_splicing_export_counts.sas

#PRG0=analysis_splicing_add_annotations.sas
# 2 hours?

## Data prep programs
#PRG0=analysis_eqtl_stack_expression.sas
	# 3 hours
#PRG0=analysis_eqtl_merge_snp_expression.sas
	# 11 hours, 1-26 PBS array
#PRG0=analysis_eqtl_prep_for_anova.sas
	# 30 min, 1-26 PBS array
#PRG0=analysis_eqtl_check_counts.sas
	# 30 min, 1-26 PBS array

## ANOVAs
#PRG0=analysis_eqtl_run_anova.sas
	# 2 hours, 1-26 PBS array
#PRG0=analysis_eqtl_export_results.sas
	# 2 hours, no array

PRG0=analysis_eqtl_calc_mean_per_genotype.sas

TMPDIR=/scratch/lfs/patcon/jnewman/sas_temp/eqtls_${PBS_ARRAYID}

if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi
cd $SASDATA

sas -log $LOGS/${PRG0}_${PBS_ARRAYID}.log -work $TMPDIR -sysin $SCRIPTS/$PRG0 -sysparm ${PBS_ARRAYID}

rm $TMPDIR/*

