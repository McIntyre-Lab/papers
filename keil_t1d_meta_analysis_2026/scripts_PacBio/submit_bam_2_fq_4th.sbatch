#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=prepMinimap.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/prepMinimap.%A_%a.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=3gb
#SBATCH --array=2-13

## 4th run!!
## convert bam to fq for input into minimap2

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio

RUN=fourth
RUNID=m64047_210312_005125

IND=$PROJ/isoseq3_refine/${RUN}_run

OUTD=${IND}/bam2fq
    mkdir -p $OUTD

module load samtools/1.10
module load bedtools/2.29.2

## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_barcode_2nd_3rd.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
BC=${ARRAY[2]}

bedtools bamtofastq \
    -i $IND/${SAMPLEID}_${RUNID}_demuxed.${BC}.flnc.bam \
    -fq $OUTD/${SAMPLEID}_${RUNID}_demuxed.${BC}.flnc.fastq
