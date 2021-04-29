#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
TAPPAS=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter

## For genes that are DE in at least one genotype, determine
##     how many are consistently up or down in FC
##     (>0, regarless of significance)

python $PROJ/scripts/count_DE_FC_consistency.py \
    -i ${TAPPAS}/all_genotype_flag_tappas_DE_gene_results_w_GO.csv \
    -g ${TAPPAS}/all_genotype_DE_FC_consistency.csv \
    -o ${TAPPAS}/all_genotype_flag_tappas_DE_consistency_w_GO.csv \
    > ${TAPPAS}/all_genotype_DE_FC_consistency_counts.txt
