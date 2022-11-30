#!/bin/sh

## Plot bar charts of sex biased expression
## For all expressed genes in each species on X and autosomes

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
Rscript ${SCRIPTS}/chip_rna_X_A_barcharts_expression_trend_only.R \
    ${MEL} ${SIM} ${ORTHO} ${OUTD}/mel_sim_X_A_gene_trend mel sim \
    > ${OUTD}/X_A_gene_trend_plot_counts.txt
