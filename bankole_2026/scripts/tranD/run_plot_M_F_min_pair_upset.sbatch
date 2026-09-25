#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=plot_M_F
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/plot_M_F_%A_%a.out
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --array=3,6

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Get directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
IND=${PROJ}/trand_sexBias



## Get design file information
DESIGN_FILE=$PROJ/design_files/design_species_tag_directory.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
TAG=${ARRAY[1]}


OUTD=${IND}/tranD_d${SPECIES}_M_vs_F
IN_FILE=${OUTD}/d${SPECIES}_male_vs_d${SPECIES}_female_minimum_pairwise_transcript_distance.csv
ROZ=${IND}/roz_d${SPECIES}
    mkdir -p ${ROZ}


## Make min pair upset plot of male biased and female biasedtranscripts


## Subset TranD pairwise output for only rows where:
##  [ flag_min_match_dmel_male OR flag_min_match_dmel_female]
awk -F',' 'NR==1 || $68==1 || $70==1' ${IN_FILE} > ${ROZ}/d${SPECIES}_M_F_minimum_pairwise_distance.csv


## Make smush upset plot of minimum pairs
python $PROJ/scripts/ignore_AS_type_02ndk.py \
    -i ${ROZ}/d${SPECIES}_M_F_minimum_pairwise_distance.csv \
    -i3 -i5 \
    -o ${OUTD}
