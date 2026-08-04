#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=espK
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/espK_%A-%a.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=126gb
#SBATCH --account=mcintyre
#SBATCH --array=1

module purge

## Description

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

PROJ=/blue/mcintyre/share/kopp_lmm_head_data

OUTD=/blue/mcintyre/share/sex_specific_splicing/axk_esp_output
	mkdir -p $OUTD

DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo $SPECIES
echo $GENOME

ES_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_${GENOME}_ujc_es.gtf
ER_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_${GENOME}_ujc_er.gtf

python $UTIL/compare_gtf_to_ES_01ksb.py \
	-es $ES_GTF \
	-er $ER_GTF \
	-i  $PROJ/${SPECIES}_data_2_${GENOME}_ujc_noMultiGene.gtf \
	-o  $OUTD
