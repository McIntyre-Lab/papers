#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=koppCatXLOC
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/catXLOC_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Description: cat TR XLOC output

## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
DIR=$PROJ/gffcompare_by_TR

echo "sampleID,numRead,numXLOCRead,numNonXLOCRead,pctXLOCRead,numGene,numXLOCGene,numNonXLOCGene" > $DIR/dser_data_2_dser1_bin_XLOC_cnts.csv
cat $DIR/dser_*_rep*_TR*_XLOC_cnts.csv >> $DIR/dser_data_2_dser1_bin_XLOC_cnts.csv
rm $DIR/dser_*_rep*_TR*_XLOC_cnts.csv
