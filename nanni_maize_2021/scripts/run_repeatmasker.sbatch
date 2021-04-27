#!/bin/bash
#SBATCH --job-name=repeatMasker
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --qos=mcintyre
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 24:00:00
#SBATCH -o /ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB/scripts/SLURM_LOGS/repeatMasker2_%A.out

module purge
module load repeatmasker

## Annotate repeat regions

## Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
PACBIO=/ufrc/mcintyre/share/maize_ainsworth/sqanti_post_filter
OUTPUT=$PROJ/repeatmasker
    mkdir -p ${OUTPUT}

## Get species name and fasta file for RepeatMasker
NAME=maize
SPECIES=maize
FA=${PACBIO}/sqanti_filtered_corrected.fasta

## cd to output directory for temporary files to be placed there
cd ${OUTPUT}

## Run repeatmasker
RepeatMasker ${FA} -s -species ${NAME} -html -gff -dir ${OUTPUT} -pa 8
