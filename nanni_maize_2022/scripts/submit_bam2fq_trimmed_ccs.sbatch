#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=bam2fq
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/bam2fq_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-11
#SBATCH --dependency=afterok:8311717

## convert bam after lima to fastq files
##     for aligning to reference genome with minimap

module load samtools/1.10 bedtools/2.29.2

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/check_isoseq3_lost_genes/bam2fq
    mkdir -p ${OUTD}
IND=$PROJ/isoseq3_analysis/library_processing

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

## Set input directory to the genotype of the sample
INGENO=${IND}/${GENO}

date

echo "
*** SAM to FQ***

Input: ${INGENO}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.bam
"
bedtools bamtofastq \
    -i ${INGENO}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.bam \
    -fq ${OUTD}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.fq

date
