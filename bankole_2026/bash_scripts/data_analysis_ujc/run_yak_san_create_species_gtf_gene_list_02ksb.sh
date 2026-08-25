#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rlrGnLst
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/gnLst_%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=1-2

PROJ=/blue/mcintyre/share/transcript_ortholog

## Get design file information (1-2)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

DESIGN_FILE=$PROJ/design_files
OUT=$PROJ/gffcompare_by_TR

awk '{print $10}' $OUT/${SPECIES}_data_2_${GENOME}_catTR.gtf | sort | uniq | tr -d '";' > $DESIGN_FILE/list_${SPECIES}_data_gene.csv

wc -l $DESIGN_FILE/list_${SPECIES}_data_gene.csv
