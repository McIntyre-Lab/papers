#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=5Sp_esp_full
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/5Sp_esp_full_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --account=mcintyre
#SBATCH --array=0-4

module purge

## Description

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

PROJ=/blue/mcintyre/share/sex_specific_splicing/

IND=$PROJ/fiveSpecies_annotations
OUTD=$IND

GENOME_ARR=(dmel6 dsim2 dsan1 dyak2 dser1)

GENOME=${GENOME_ARR[${SLURM_ARRAY_TASK_ID}]}

ANNO_FILE=$IND/fiveSpecies_2_${GENOME}_ujc.gtf

python $UTIL/compare_gtf_to_ES_01ksb.py \
	-es $IND/fiveSpecies_2_${GENOME}_ujc_es.gtf \
	-er $IND/fiveSpecies_2_${GENOME}_ujc_er.gtf \
	-i  $ANNO_FILE \
	-o  $OUTD
