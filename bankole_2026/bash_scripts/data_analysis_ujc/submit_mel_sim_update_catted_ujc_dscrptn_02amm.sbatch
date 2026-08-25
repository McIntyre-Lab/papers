#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=ujc_dsc_RMG
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/ujc_dsc_R_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account=mcintyre
#SBATCH --array=1-2

module purge

## Description: Fix UJC dscrptn file post-catting. add flagMultiGene
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

OUTD=$PROJ

## Get design file information (2 rows, dsan and dyak)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=d${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

python $UTIL/find_multiGene_jxnHash_in_catted_ujc_dsc_file_01ksb.py \
	-i $PROJ/${SPECIES}_data_2_${GENOME}_ujc_dscrptn_simpleCat.csv \
	-o $PROJ/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv

