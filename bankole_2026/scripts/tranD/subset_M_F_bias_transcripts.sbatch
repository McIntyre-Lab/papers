#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=subset
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/subset_%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=2,3,6

## input dir
PROJ=/blue/mcintyre/share/sex_specific_splicing
IND=${PROJ}/trand_sexBias

module load python/3.12

## create output dirs

OUTD=${PROJ}/trand_sexBias
    mkdir -p $OUTD

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_species_tag_directory.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
TAG=${ARRAY[1]}

ROZ=$OUTD/roz_${SPECIES}
    mkdir -p ${ROZ}

GTF=${PROJ}/d${SPECIES}_data_2_${TAG}_ujc_noMultiGene.gtf

##Selecting Male biased ujcs
     awk -F',' 'NR > 1 {print $1}' ${IND}/d${SPECIES}_male_ujc_4_TranD.csv \
     > ${ROZ}/list_d${SPECIES}_male_ujcs_4_TranD.txt

    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${GTF} \
        -t transcript_id \
        -i ${ROZ}/list_d${SPECIES}_male_ujcs_4_TranD.txt \
        -o ${OUTD}/d${SPECIES}_male_ujcs_4_TranD.gtf

##Selecting female biased ujcs
     awk -F',' 'NR > 1 {print $1}' ${IND}/d${SPECIES}_female_ujc_4_TranD.csv \
     > ${ROZ}/list_d${SPECIES}_female_ujcs_4_TranD.txt

    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${GTF} \
        -t transcript_id \
        -i ${ROZ}/list_d${SPECIES}_female_ujcs_4_TranD.txt \
        -o ${OUTD}/d${SPECIES}_female_ujcs_4_TranD.gtf
