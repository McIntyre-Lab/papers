#!/bin/bash

## Make sure to activate conda environment with upsetplot package
##     AVN used custom env "mylab" on grimshawi

## Plot boxplot of gene mean expression TPM values (sum of transcript means)
##     within the groups of genes
##         1) shared by Wang 2018 leaf transcriptome and
##            5 genotype/2 treatment leaf transcriptome, and
##         2) only in 5 genotype/2 treatment leaf transcriptome

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
OUTD=$PROJ/sqanti_classification_category_subset/plot_consol_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
EA1=$PROJ/sqanti_classification_category_subset/EA_annotations
EA2=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/maize_leaf_150bp_annotations

python $PROJ/scripts/plot_compare_leaf_transcriptome_expression.py \
    -f ${IND}/flag_on_off_rsem_expression_isoforms_TPM_5_full_table.tsv \
    -a ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_id_2_gene_id.csv \
    -1 ${EA1}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv \
    -2 ${EA2}/maize_leaf_consol_event2transcript2gene_index.csv \
    --missing1 $PROJ/sqanti_classification_category_subset/plot_FSM_ISM_zmtr_NIC_NNC_monoexon_filter/diff_input_ea_output_transcript.csv \
    --missing2 $PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome/maize_leaf_diff_input_ea_output_transcript.csv \
    -p ${OUTD}/compare_leaf_transcriptomes
