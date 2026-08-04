#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=ujcGtfK
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/ujcGtfK_%A-%a.out
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=56gb
#SBATCH --account=mcintyre
#SBATCH --array=1

module purge

## Description: Subset ujc GTF for nonMultiGene jxnHash.

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

## Get design file information
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

# Create list of jxnHash that are not multiGene
tail -n +2 $PROJ/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_noMultiGene.csv | awk -F, '{print $4}' | sort | uniq > $SLURM_TMPDIR/list_${SPECIES}_data_2_${GENOME}_noMultiGene_jxnHash.csv
NMGLIST=$SLURM_TMPDIR/list_${SPECIES}_data_2_${GENOME}_noMultiGene_jxnHash.csv

# Subset GTF for only multiGene jxnHash for the GTF
grep -f $NMGLIST $PROJ/${SPECIES}_data_2_${GENOME}_ujc_simpleCat.gtf > $PROJ/${SPECIES}_data_2_${GENOME}_ujc_noMultiGene.gtf
