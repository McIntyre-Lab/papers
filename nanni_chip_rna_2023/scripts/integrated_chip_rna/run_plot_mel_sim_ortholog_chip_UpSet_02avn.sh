#!/bin/sh

## Make sure to activate conda environment with upsetplot package
##     AVN used custom env "mylab" on grimshawi

export PATH=$HOME/conda/.envs/mylab/bin:$PATH

## Plot UpSet plots (X and Autosomes separately) of the
##     number/proportion of genes that have H3K4me3/H3K27me2me3
##     marks present in mel and sim

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/manuscript/supp_files
OUTPUT=${PROJ}/manuscript/figures/ortho_plots
    mkdir -p ${OUTPUT}

## Plot
python ${SCRIPTS}/plot_mel_sim_ortholog_chip_UpSet_02avn.py \
    -i ${INPUT}/dmel_dsim_ortholog_chip_rna_flags.csv \
    -l ${OUTPUT}/mel_sim_ortho_gene_chip_UpSet_plot.log \
    -p ${OUTPUT}/mel_sim_ortho_gene
