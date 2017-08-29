#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=run_sas
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/concannon/share/jnewman/scripts/PBS_LOGS/out.aln_splicing.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1

## 25 for splicing, 639 for eQTL
module load sas/9.4

#Set up variables

PROJ=/ufrc/concannon/share/jnewman
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
#PRG0=analysis_eqtl_merge_snp_expression_jrbn2.sas
	# 11 hours, 1-26 PBS array
#PRG0=analysis_eqtl_prep_for_anova.sas
	# 30 min, 1-26 PBS array
#PRG0=analysis_eqtl_check_counts.sas
	# 30 min, 1-26 PBS array

## ANOVAs
#PRG0=analysis_eqtl_run_anova.sas
	# 2 hours, 1-26 PBS array
PRG0=analysis_eqtl_export_results.sas
	# 2 hours, no array

#PRG0=analysis_eqtl_calc_mean_per_genotype_${PBS_ARRAYID}.sas

TMPDIR=/ufrc/concannon/share/jnewman/sas_temp/eqtls_${SLURM_ARRAY_TASK_ID}

if [ ! -e $TMPDIR ]; then mkdir -p $TMPDIR; fi
cd $SASDATA

sas -log $LOGS/${PRG0}_${SLURM_ARRAY_TASK_ID}.log -work $TMPDIR -sysin $SCRIPTS/$PRG0 -sysparm ${SLURM_ARRAY_TASK_ID}

rm $TMPDIR/*

