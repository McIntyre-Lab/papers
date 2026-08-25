#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rlrCatGTF
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/catGTF_%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --array=1-2

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: cat noXLOC gtf into one big species GTF

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
DIR=$PROJ/gffcompare_by_TR

## Get design file information (1-2)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

cat $DIR/${SPECIES}_*_TR*_2_${GENOME}_noXLOC.gtf > $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_catTR.gtf

# Sort by seqname, start, and transcriptID
# sort -t $'\t' -k1,1 -k4,4n -k10,10 $SLURM_TMPDIR/${SPECIES}_data_2_${GENOME}_catTR.gtf > $DIR/${SPECIES}_data_2_${GENOME}_catTR.gtf
