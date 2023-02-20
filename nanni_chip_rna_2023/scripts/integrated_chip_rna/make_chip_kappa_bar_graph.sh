#!/bin/sh

## Plot Cohen's kappa between males and females
##     in H3K4me3 and H3K27me2me3 of each feature

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/ChIPseq/detection_above_background/features
OUTPUT=${PROJ}/figures/kappa_bar_graph
    mkdir -p ${OUTPUT}

## Get input files
MEL=${INPUT}/mel_chip_features_flag_kappas.csv
SIM=${INPUT}/sim_chip_features_flag_kappas.csv

## Plot
Rscript ${SCRIPTS}/chip_kappa_bar_graph.R \
    ${MEL} ${SIM} ${OUTPUT}/mel_sim_kappa_bar_graph \
    mel sim
