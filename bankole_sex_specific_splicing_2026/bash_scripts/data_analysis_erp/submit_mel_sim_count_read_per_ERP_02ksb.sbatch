#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=erpCntR
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/erpCntR_%A-%a.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --account=mcintyre
#SBATCH --array=1-2

module purge

## Description

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

PROJ=/blue/mcintyre/share/rmg_lmm_dros_data

DIR=/blue/mcintyre/share/sex_specific_splicing/rmg_erp_output

DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=d${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

python $UTIL/count_reads_per_ERP_01ksb.py \
	-i $DIR/fiveSpecies_2_${GENOME}_ujc_er_vs_${SPECIES}_data_2_${GENOME}_ujc_noMultiGene_infoERP.csv \
	-c $PROJ/${SPECIES}_data_2_${GENOME}_ujc_count.csv \
	-o $DIR \
	-x fiveSpecies_2_${GENOME}_ujc_er_vs_${SPECIES}_data_2_${GENOME}_ujc_noMultiGene
