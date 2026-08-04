#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=espCntK
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/espCntK_%A-%a.out
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=26gb
#SBATCH --account=mcintyre
#SBATCH --array=1
#SBATCH --dependency=afterok:49310073

module purge

## Description

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

PROJ=/blue/mcintyre/share/kopp_lmm_head_data

DIR=/blue/mcintyre/share/sex_specific_splicing/axk_esp_output

DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

python $UTIL/count_reads_per_ESP_01ksb.py \
	-i $DIR/fiveSpecies_2_${GENOME}_ujc_es_vs_${SPECIES}_data_2_${GENOME}_ujc_noMultiGene_infoESP.csv \
	-c $PROJ/${SPECIES}_data_2_${GENOME}_ujc_count.csv \
	-o $DIR \
	-x fiveSpecies_2_${GENOME}_ujc_es_vs_${SPECIES}_data_2_${GENOME}_ujc_noMultiGene
