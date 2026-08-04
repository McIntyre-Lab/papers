#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=trand_gene_5s
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/trand_gene_5s_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=0-4
#SBATCH --qos=mcintyre-b

module purge

## Description

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/

IND=$PROJ/fiveSpecies_annotations

ARRAY=(dmel6 dsim2 dsan1 dyak2 dser1)

GENOME=${ARRAY[$SLURM_ARRAY_TASK_ID]}

OUTD=$PROJ/trand_1gtf_geneMode_fiveSpecies
	mkdir -p $OUTD

date

NAME=fiveSpecies_2_${GENOME}_ujc

trand \
$IND/${NAME}.gtf \
-o $OUTD \
-e gene \
-x $NAME \
-f \
-n $SLURM_CPUS_ON_NODE \
--verbose \
-l $OUTD/fiveSpecies_2_${GENOME}_ujc.log

date
