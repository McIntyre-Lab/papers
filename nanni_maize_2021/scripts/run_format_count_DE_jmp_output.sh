#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge tappAS output files with detection flags

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter

for TYPE in biological_process cellular_component molecular_function; do
    python $PROJ/scripts/format_count_DE_jmp_output.py \
        --JMPfile ${IND}/GO_gse_tappas_DE_${TYPE}.csv \
        --GOterm $PROJ/go_enrichment/GO_table.txt \
        -o ${IND}/GO_gse_tappas_DE_${TYPE} \
        > ${IND}/GO_gse_tappas_DE_${TYPE}_counts.txt
done
