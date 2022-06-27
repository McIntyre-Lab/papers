#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_ttest_flags_with_anno.csv
SIM=${INPUT}/sim_ttest_flags_with_anno.csv

for SPECIES in mel sim; do

    if [[ ${SPECIES} == "mel" ]]; then
        ANNOT=${MEL}
    else
        ANNOT=${SIM}
    fi

    ## Get frequencies of rna detected APN > 5 with gene_sex_bias
    ## gene_sex_bias is based off of APN > 0 detection with 2-fold change
    python3 ${SCRIPTS}/count_gene_sex_bias_ttest_foldchange_rna_03avn.py \
        -i ${ANNOT} \
        -c ${INPUT}/${SPECIES}_gene_sex_bias_counts.txt \
        -o ${INPUT}/${SPECIES}_ttest_flags_with_anno_combo_flag.csv
done
