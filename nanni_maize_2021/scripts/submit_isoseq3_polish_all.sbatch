#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_polish
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/pol_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-264


### isoseq3 analysis of maize test PacBio data - polish by individual

### note script uses 2 design files!

## 19 is array design 49-72
## 113 is array design 1-24 (c123 amb) 
## 46 is array design 145-168
## 89 is array design 217-240

## 21 is 97-120
## 21-2 is 73-96
## 120 is 26-48
## 42 is 121-144
## 67 is 169-192
## 70 is 193-216
## 96 is 241-264

module load isoseq3/3.0.0 pacbio/5.1.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/isoseq3_analysis/polish_by_individual
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IND=$PROJ/isoseq3_analysis/cluster_by_individual

## Get info from first design file

DESIGN_FILE=$PROJ/design_files/design_maize_PacBio_cluster_test_unpolished_path_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

ID=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
UNPOLISHED=${ARRAY[3]}

## Set input and output directories to the correct individual sample
OUT=$OUTD/$GENO/$TRT
    mkdir -p $OUT
IN=$IND/$GENO/$TRT

## Get subreads.bam input file from a different design file
## identify subreads.bam for ID, GENO, TRT combinations
## note that df below is sampleID specific!
SUBREAD=$(ls $(grep -e $ID -e $GENO -e $TRT $PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_${ID}.csv | awk 'BEGIN{FS=","}{print $5}')/*.subreads.bam)
echo "Subread Input File : $SUBREAD
"

## Get bam file basename and number
SAMPLE_NAME=$(echo $UNPOLISHED | awk 'BEGIN{FS="."}{print $1}')
NUM=$(echo "$UNPOLISHED" | awk -F "." '{print $(NF-1)}')
date
echo "***Polish***

Inputs : 
$IN/$UNPOLISHED 
$SUBREAD
"
## Polish
## creates a consensus sequence for each clustered transcripts
## input is a unpolished.bam or a split (from isoseq3 cluster --split-bam) unpolished.#.bam
## outputs polished.transcriptset.#.bam

#isoseq3 polish $IN/$UNPOLISHED $SUBREAD $OUT/${SAMPLE_NAME}.polished.transcriptset.${NUM}.bam

