#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=koppGnLst
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/gnLst_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb

PROJ=/blue/mcintyre/share/kopp_lmm_head_data

DESIGN_FILE=$PROJ/design_files
OUT=$PROJ/gffcompare_by_TR

awk '{print $10}' $OUT/dser_data_2_dser1_catTR.gtf | sort | uniq | tr -d '";' > $DESIGN_FILE/list_dser_data_gene.csv

wc -l $DESIGN_FILE/list_*_data_gene.csv
