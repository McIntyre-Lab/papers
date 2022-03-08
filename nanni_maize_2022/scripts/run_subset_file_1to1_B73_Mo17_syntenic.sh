#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
INB=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/tappas_output
SYN=~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data/synteny_Mo17Cau_B73v4

## Subset tappAS output merged with GO for B73-Mo17 1-to-1 syntenic pairs

python $PROJ/scripts/subset_file_1to1_B73_Mo17_syntenic.py \
    -i ${INB}/all_gene_flag_tappas_DE_gene_results_w_GO.csv \
    -s ${SYN}/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv \
    -g B73 \
    -c gene_id \
    -o ${INB}/all_gene_flag_tappas_DE_gene_results_w_GO_1to1BMo.csv

