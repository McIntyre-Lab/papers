#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_lib_proc
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/iso_lib_proc_%A_%a.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11

### isoseq3 analysis of maize samples- library processing
###
### starting input is a *.subreads.bam file - sequences include subreads separated by primer and adapter sequences

### modify DF to do different sampleSets


module load isoseq3/3.0.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/isoseq3_analysis/library_processing
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IN=$PROJ/original_data
ROZ=$PROJ/roz_isoseq3
    if [ ! -e $ROZ ]; then mkdir -p $ROZ; fi

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SUBREADS_BAM=$IN/${SAMPLE_PATH}/*.subreads.bam
PRIMERS=$PROJ/pacbio_primers.fasta
SAMPLE_NAME=${ID}_${GENO}_${TRT}

## Set output directory to the species of the sample
OUT=$OUTD/$GENO
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
LOG=$OUTD/log
    if [ ! -e $LOG ]; then mkdir -p $LOG; fi

date

echo "
***Consensus Calling***

Input: ${SUBREADS_BAM}
"
## Consensus Calling
## outputs consensus sequence from each ZMW with at least one full pass in ccs.bam
## --noPolish only output original template derived from the partial order alignment of subreads
## --minPasses minimum number of full passes per read (each primer is seen at least once)

ccs --force --noPolish --minPasses 1 --numThreads 12 --logFile=$LOG/${SAMPLE_NAME}.log ${SUBREADS_BAM} $OUT/${SAMPLE_NAME}.ccs.bam

echo "Output: $OUT/${SAMPLE_NAME}.ccs.bam
"

date

echo "
***Primer Removal***

Input: $OUT/${SAMPLE_NAME}.ccs.bam
"
## Primer Removal (and Demultiplex)
## trim primers (and barcodes if used)
## --isoseq IsoSeq mode
## --no-pbi do not generate .bam.pbi output

lima --isoseq --no-pbi --num-threads 12 $OUT/${SAMPLE_NAME}.ccs.bam $PRIMERS $OUT/${SAMPLE_NAME}.ccs.trimmed.bam

rm -rf $ROZ
