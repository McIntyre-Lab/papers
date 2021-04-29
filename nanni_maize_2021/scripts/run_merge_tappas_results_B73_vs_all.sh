#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge tappAS output files with detection flags

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
TAPPAS=${IND}/tappas_output_TMM_norm_1CPMfilter/B73_Amb_vs_all_Amb
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info

python $PROJ/scripts/merge_tappas_B73_Amb_vs_all_DIU_results.py \
    -t ${TAPPAS} \
    -f ${IND}/flag_on_off_rsem_expression_isoforms_TPM_5.tsv \
    -a ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_id_2_gene_id.csv \
    -e B73_P1_C7_Ele -e B73_P4_C1_Amb -e B73_P4_C6_Amb \
    -o ${TAPPAS}/B73_vs_all_flag_tappas_DIU_gene_results.csv \
    -m ${TAPPAS}/B73_vs_all_atLeast1geno_DIU_majorIsoformSwitch.csv \
    > ${TAPPAS}/B73_vs_all_DIU_gene_counts.txt
