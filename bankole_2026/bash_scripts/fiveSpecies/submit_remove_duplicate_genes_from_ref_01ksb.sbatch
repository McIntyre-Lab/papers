#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=dupeGene
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/dupeGene_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --account=mcintyre
#SBATCH --array=2-6

#module purge

## Description: Find duplicate genes from a reference annotation and remove them (if they exist)

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing/

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

## Get design file information
DESIGN_FILE=$PROJ/design_files/df_run_dupe_gene_01ksb.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

ANNONAME=${ARRAY[0]}
REFFOLDER=${ARRAY[1]}
ANNO=${ARRAY[2]}

python $UTIL/find_duplicated_genes_in_gtf_01ksb.py \
	-a $ANNO \
	-p $ANNONAME \
	-o $REFFOLDER

