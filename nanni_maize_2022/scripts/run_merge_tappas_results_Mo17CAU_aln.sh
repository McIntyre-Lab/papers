#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge tappAS output files with detection flags

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/tappas_output_mo17_cau_cvr_cnt
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info

## Merge tappas results from each genotype and make flags
python $PROJ/scripts/merge_flag_tappas_DE.py \
    -d ${IND} \
    -o ${IND}/all_gene_flag_tappas_DE_gene_results.csv \
    > ${IND}/all_gene_flag_tappas_DE_gene_result.log

## Merge with GO terms
python $PROJ/scripts/merge_GO_tappas_results.py \
    -t ${IND}/all_gene_flag_tappas_DE_gene_results.csv \
    --type DE \
    --GOid ${REF}/RefGenV4/sasdata/GO_IDs_catted.csv \
    --GOterm $PROJ/go_enrichment/GO_table.txt \
    --output-full ${IND}/all_gene_flag_tappas_DE_gene_results_w_GO.csv \
    --output-5 ${IND}/list_DE_all_5_gene_w_GO.csv
