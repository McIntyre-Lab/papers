#!/bin/sh


#########################
###   Peak Analysis   ###
#########################

## Call peaks using MACS2 and evaluate replicates


## Set Output directories
export MACS=$PROJ/chipseq/macs2_peak_calling
    mkdir -p ${MACS}
export COMBINE=$PROJ/chipseq/macs2_consensus_peaks
    mkdir -p ${COMBINE}
export ANNOTATE=$PROJ/chipseq/annotate_macs2_consensus
    mkdir -p ${ANNOTATE}
export PILES=${ALN}/mpileups
    mkdir -p ${PILES}
export CVR=$PROJ/chipseq/coverage_cnts/macs2_consensus_peaks
    mkdir -p ${CVR}/logs
export DABG=$PROJ/chipseq/detection_above_background/macs2_consensus_peaks
    mkdir -p ${DABG}
export FRIP=$PROJ/chipseq/frip/macs2_consensus_peaks
    mkdir -p ${FRIP}

## Get BAM file directory
export BAMFILES=${ALN}/bam_files


## Call MACS2 peaks for each sample
#source ${SCRIPTS}/call_MACS2_peaks.sh


## Combine replicate peaks
#source ${SCRIPTS}/combine_rep_MACS2_peaks.sh


## Annotate MACS2 consensus peaks
#source ${SCRIPTS}/annotate_peaks_03avn.sh


## Coverage counts on combined MACS2 peaks
#source ${SCRIPTS}/coverage_counts_peaks.sh


## Detection above background (DABG) of MACS2 peaks
#source ${SCRIPTS}/detection_above_background_peaks.sh


## Fraction of reads in peaks (FRiP)of MACS2 peaks
#source ${SCRIPTS}/frip_peaks.sh

## Plot FRiP values in box and whisker plot, split by ChIP marks
module load R/3.6
SUMMARY=${FRIP}/summary_frip_replicates.csv
BOXPLOT=${FRIP}/summary_frip_replicates
Rscript ${SCRIPTS}/boxplot_frip_peaks.R \
    ${SUMMARY} \
    ${BOXPLOT} \
    ${SPECIES1} ${SPECIES2}
