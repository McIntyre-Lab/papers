#!/bin/sh

## Plot horizontal stacked bar charts of the proportion of genes
##     on the D. melanogaster and D. simulans X and Autosomes
##     based on the presence of H3K4me3 or H3K27me2me3

## Include plot of subset to one-to-one orthologs within each species

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
SCRIPTS=$PROJ/scripts/integrated_chip_rna
IND=$PROJ/manuscript/supp_files
OUTD=$PROJ/manuscript/figures/X_A_plots/new_plots
    mkdir -p ${OUTD}

## Get input files from supplementary
MEL=${IND}/dmel_chip_rna_flags.csv
SIM=${IND}/dsim_chip_rna_flags.csv
ORTHO=${IND}/dmel_dsim_ortholog_chip_rna_flags.csv

## Plot gene-level
Rscript ${SCRIPTS}/chip_presence_X_A_barcharts_03avn.R \
    ${MEL} ${SIM} ${ORTHO} ${OUTD}/mel_sim_X_A_chip_presence mel sim \
    > ${OUTD}/mel_sim_X_A_chip_presence_counts.txt
