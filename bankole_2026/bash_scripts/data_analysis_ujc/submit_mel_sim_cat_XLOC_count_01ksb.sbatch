#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rmgCatXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/catXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=1-2

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: cat TR XLOC output

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
DIR=$PROJ/gffcompare_by_TR

## Get design file information (1-2)
DESIGN_FILE=$PROJ/design_files/list_species_geno.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
GENOME=${ARRAY[1]}

echo "sampleID,numRead,numXLOCRead,numNonXLOCRead,pctXLOCRead,numGene,numXLOCGene,numNonXLOCGene" > $DIR/${SPECIES}_data_2_${GENOME}_bin_XLOC_cnts.csv
cat $DIR/${SPECIES}_*_rep*_TR*_XLOC_cnts.csv >> $DIR/${SPECIES}_data_2_${GENOME}_bin_XLOC_cnts.csv
rm $DIR/${SPECIES}_*_rep*_TR*_XLOC_cnts.csv
