#!/bin/sh


export PATH=$HOME/conda/.envs/mylab/bin:$PATH
#export PATH=/Users/adalena/opt/anaconda3/envs/mylab/bin:$PATH

## Make barplots of 4 panels (male-biased, female-biased, and X, A) 
##     with 8 bars each (open chromatin: mel presence vs. not,
##     sim presence vs. not; closed chroamtin: mel presence vs. not,
##     sim presence vs. not)

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/supp_files
OUTPUT=${PROJ}/figures/exp_chrom_barplots
#PROJ=~/Desktop/dros_figures
#SCRIPTS=$PROJ/scripts/integrated_chip_rna
#INPUT=$PROJ/supp_files
#OUTPUT=$PROJ/figures/exp_chrom_barplots
    mkdir -p ${OUTPUT}

## Plot
python ${SCRIPTS}/exp_chrom_barplot_XA_02avn.py \
    -m ${INPUT}/dmel_chip_rna_flags.csv \
    -s ${INPUT}/dsim_chip_rna_flags.csv \
    -o ${INPUT}/dmel_dsim_ortholog_chip_rna_flags.csv \
    -d ${OUTPUT}
