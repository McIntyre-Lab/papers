#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge tappAS output files with detection flags

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/tappas_output_mo17_cau_cvr_cnt

## Merge with GO terms
for TYPE in biological_process cellular_component molecular_function; do
    python $PROJ/scripts/format_count_DE_jmp_output.py \
        --GOterm $PROJ/go_enrichment/GO_table.txt \
        --JMPfile ${IND}/GO_gse_tappas_DE_${TYPE}.csv \
        -o ${IND}/GO_gse_tappas_DE_${TYPE} \
        > ${IND}/GO_gse_tappas_DE_${TYPE}_counts.txt
done
