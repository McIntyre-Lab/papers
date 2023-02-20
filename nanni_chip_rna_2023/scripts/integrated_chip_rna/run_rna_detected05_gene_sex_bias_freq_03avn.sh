#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_ttest_flags_with_anno_combo_flag.csv
SIM=${INPUT}/sim_ttest_flags_with_anno_combo_flag.csv

for SPECIES in mel sim; do

    if [[ ${SPECIES} == "mel" ]]; then
        ANNOT=${MEL}
    else
        ANNOT=${SIM}
    fi

    ## Get frequencies of rna detected APN > 5 with gene_sex_bias
    ## gene_sex_bias is based off of APN > 0 detection with 2-fold change
    ## and split by sex-limited and sex-biased genes
    python ${SCRIPTS}/rna_detected05_gene_sex_bias_freq_05avn.py \
        -i ${ANNOT} \
        -o ${INPUT}/${SPECIES}_gene_sex_bias_detected05_counts.txt

done
