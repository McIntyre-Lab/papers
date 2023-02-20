#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=annotation
#SBATCH --output=/ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/annotation_pipeline_%A.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

### Set Global Variables

## Project and scripts directories
export PROJ=/ufrc/mcintyre/share/etoh_srna
export SCRIPTS=$PROJ/scripts/annotation

## Species Names (will be used in output file names)
export SPECIES1="mel"
export SPECIES2="sim"

## Corresponding FASTA and GTF files for each species
export FASTA1=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
export GTF1=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-r6.17.gtf
export FASTA2=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
export GTF2=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-r2.02.gtf

## Corresponding Event Analysis annotation prefix
export EA1=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_annotations_100bp_reads/deml617_100bp
export EA2=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_100bp_reads/dsim202_100bp

###############################
### Collect unique features ###
###############################
source ${SCRIPTS}/unique_features_04avn.sh
