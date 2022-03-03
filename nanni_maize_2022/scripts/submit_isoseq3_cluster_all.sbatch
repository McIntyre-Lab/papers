#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_clus_indiv
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/iso_clus_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11

### isoseq3 analysis of maize test PacBio data - cluster by individual 
### use polyA options

module load isoseq3/3.0.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/isoseq3_analysis/cluster_by_individual
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IND=$PROJ/isoseq3_analysis/library_processing

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

## Set input and output directories to the correct species
OUT=$OUTD/$GENO/$TRT
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IN=$IND/$GENO

SAMPLE_NAME=${ID}_${GENO}_${TRT}

date

echo "
***Cluster By Individual***
"

## Cluster
## input is the consensus read bam *.ccs.bam
## removes polyA and artificial concatemers to create FLNC (Full-Length,Non-Concatemer) reads, then puts them in transcript sets
## --split-bam 24 splits the output into 24 files for downstream parallelization
## outputs unpolished.#.bam and unpolished.transcriptset.#.xml where # is the file number after splitting
## use --require-polya --> trims polyA seq

isoseq3 cluster --require-polya --num-threads 12 -v --split-bam 24 $IN/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.bam $OUT/${SAMPLE_NAME}.unpolished.bam
