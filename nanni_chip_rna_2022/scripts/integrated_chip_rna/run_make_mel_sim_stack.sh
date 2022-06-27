#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_ttest_flags_with_anno_combo_flag.csv
SIM=${INPUT}/sim_ttest_flags_with_anno_combo_flag.csv
ORTHO=$PROJ/CHIP_RNA_ms/results/data_files/mel_sim_ortho_combo_flags.csv

python ${SCRIPTS}/make_mel_sim_stack.py \
    -m ${MEL} \
    -s ${SIM} \
    --ortholog ${ORTHO} \
    -o ${INPUT}/mel_sim_rna_chip_flag_stack.csv
