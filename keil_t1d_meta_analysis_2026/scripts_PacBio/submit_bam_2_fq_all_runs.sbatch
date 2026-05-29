#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=bam2fq.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/bam2fq.%A_%a.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/bam2fq.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --array=2-61

## convert bam to fq for input into minimap2

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio

## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_bc_2_run.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]} 
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}
SAMPLEID=${ARRAY[3]}
BC=${ARRAY[4]}



IND=$PROJ/isoseq3_refine_by_sample_isoseq/${RUN}

OUTD=${IND}/bam2fq
    mkdir -p $OUTD

module load samtools/1.15
module load bedtools/2.29.2

bedtools bamtofastq \
    -i $IND/${SAMPLEID}_${RUN}.flnc.bam \
    -fq $OUTD/${SAMPLEID}_${RUN}.flnc.fastq

