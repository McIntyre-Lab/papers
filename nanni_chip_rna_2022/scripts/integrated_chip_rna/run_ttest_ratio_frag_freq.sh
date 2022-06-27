#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_frag_flags_kitchen_sink.csv
SIM=${INPUT}/sim_frag_flags_kitchen_sink.csv

for SPECIES in mel sim; do

    if [[ ${SPECIES} == "mel" ]]; then
        ANNOT=${MEL}
    else
        ANNOT=${SIM}
    fi

    ## Get frequencies of fragment flags of:
    ##   1) ttest significance flags with ratio != 1
    ##   2) fold change >= 2
    ##   3) ttest significant flag and fold change >= 2
    python ${SCRIPTS}/ttest_ratio_frag_freq_02avn.py \
        -i ${ANNOT} \
        -o ${INPUT}/${SPECIES}_ttest_ratio_frag_freq.txt

done
