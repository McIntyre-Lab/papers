#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_lib_proc
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/iso_lib_proc_MaizeWang2018_%A_%a.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-107

### isoseq3 analysis of MaizeWang2018 samples (pool 1) - library processing
### have 107 *bas.h5 files with 222 *.subreads.fastq files that are barcoded
### convert *.bas.h5/*.bax.h5 legacy file formats to *.subreads.bam files
### starting input for ccs is a *.subreads.bam file - sequences include subreads separated by primer and adapter sequences
### lima will remove primers and demultiplex

module load bax2bam/0.0.9 isoseq3/3.0.0

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/isoseq3_library_processing
BAMOUT=${OUTD}/subreads_bam
    mkdir -p ${BAMOUT}
CCSOUT=${OUTD}/ccs
    mkdir -p ${CCSOUT}
LIMAOUT=${OUTD}/lima
    mkdir -p ${LIMAOUT}
LOG=${OUTD}/log
    mkdir -p ${LOG}

## Get barcode fasta file
BARCODES=$PROJ/design_files/MaizeWang2018_barcode_primer.fa

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_leaf_subread_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BAS=${ARRAY[0]}
NAME=${ARRAY[1]}

## Get data directory
DATADIR=$(dirname ${BAS})

date

echo "
*** Convert *.bas.h5 to *.subreads.bam ***
"
## Convert legacy *.bas.h5/*.bax.h5 file formats to *.subreads.bam
##     for proper input format to ccs
if [[ ! -e ${BAMOUT}/${NAME}.subreads.bam ]]; then
    bax2bam --subread \
        -o ${BAMOUT}/${NAME} \
        ${DATADIR}/${NAME}*.bax.h5
else
    echo "Subreads BAM file found: skipping bax2bam step"
fi
date

echo "
*** Consensus Calling ***
"
## Consensus Calling
## outputs consensus sequence from each ZMW with at least one full pass in ccs.bam
## --noPolish only output original template derived from the partial order alignment of subreads
## --minPasses minimum number of full passes per read (each primer is seen at least once)
if [[ ! -e ${CCSOUT}/${NAME}.ccs.bam ]]; then
    ccs --force \
        --noPolish \
        --minPasses 1 \
        --numThreads 12 \
        --logFile ${LOG}/${NAME}.log \
        ${BAMOUT}/${NAME}.subreads.bam \
        ${CCSOUT}/${NAME}.ccs.bam
else
    echo "CCS BAM file found: skipping ccs step"
fi
date

echo "
*** Demultiplex and Primer Removal ***

Input: ${CCSOUT}/${NAME}.ccs.bam
"
## Primer Removal (and Demultiplex)
## trim primers (and barcodes if used)
## --isoseq activate specialized isoseq mode for primer removal and demultiplexing
## --peek-guess Try to infer the used barcodes subset, by peeking at the first 50,000 ZMWs, whitelisting barcode pairs with more than 10 counts and mean score ≥ 45.
lima --num-threads 12 \
    --isoseq \
    --no-pbi \
    ${CCSOUT}/${NAME}.ccs.bam \
    ${BARCODES} \
    ${LIMAOUT}/${NAME}_ccs_demux.bam
date
