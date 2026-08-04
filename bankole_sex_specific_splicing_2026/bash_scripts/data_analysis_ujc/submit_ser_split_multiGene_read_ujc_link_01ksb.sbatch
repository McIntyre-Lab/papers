#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=ujc_lnk_TO
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/ujc_lnk_K_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256gb
#SBATCH --account=mcintyre
#SBATCH --array=1

module purge

## Description: Fix UJC dscrptn file post-catting. add flagMultiGene
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

PROJ=/blue/mcintyre/share/kopp_lmm_head_data
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

OUTD=$PROJ

## Get design file information (2 rows, dsan and dyak)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

python $UTIL/output_multiGene_jxnHash_using_ujc_link_file_01ksb.py \
	-d $PROJ/${SPECIES}_data_2_${GENOME}_ujc_dscrptn.csv \
	-l $PROJ/${SPECIES}_data_2_${GENOME}_ujc_xscript_link_simpleCat.csv \
	-x ${SPECIES}_data_2_${GENOME} \
	-o $PROJ

