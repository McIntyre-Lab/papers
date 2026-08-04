#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=ujc_cnt_TO
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/ujc_cnt_TO_%A-%a.out
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=1-2

module purge

## Description: Count nonMultiGene read per jxnHash.

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/transcript_ortholog
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

## Get design file information
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

# Use utility to count nonMultiGene read per jxnHash
python $UTIL/create_ujc_counts_from_catted_link_file_01ksb.py \
	-l $PROJ/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_noMultiGene.csv \
	-o $PROJ \
	-x ${SPECIES}_data_2_${GENOME}

