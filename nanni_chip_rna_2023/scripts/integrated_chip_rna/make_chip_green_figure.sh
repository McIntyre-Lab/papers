#!/bin/sh

## Plot bar charts of sex biased expression and H3K4me3

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/CHIP_RNA_ms/results/data_files
OUTPUT=${PROJ}/CHIP_RNA_ms/figures/green_fig
    mkdir -p ${OUTPUT}

## Get input files
MEL=${INPUT}/mel_chip_rna_frag_flags_anno.csv
SIM=${INPUT}/sim_chip_rna_frag_flags_anno.csv

for SPECIES in mel sim; do
    if [[ ${SPECIES} == "mel" ]]; then
        ANNOT=${MEL}
    else
        ANNOT=${SIM}
    fi
    ## Get frequencies of histone detection groups
    cut -d "," -f 1,11-14 ${ANNOT} | awk 'NR!=1' | sort | uniq -c | \
        awk 'BEGIN{print "featureType,fK4,fK27,mK4,mK27,freq"}{print $2","$1}' \
        > ${OUTPUT}/${SPECIES}_green_fig_counts.csv

    ## Plot
    Rscript ${SCRIPTS}/chip_green_figure_plot.R \
        ${OUTPUT}/${SPECIES}_green_fig_counts.csv \
        ${OUTPUT}/${SPECIES}_green_fig
done
