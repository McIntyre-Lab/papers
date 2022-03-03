#!/bin/bash

## Make sure to activate conda environment with upsetplot package
##     AVN used custom env "mylab" on grimshawi

## Plot groups of DE across the 5 genotypes where genes
##     are detected in at least one treatment to be DE

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
TAPPAS=${IND}/tappas_output_TMM_norm_1CPMfilter
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info/RefGenV4

#python $PROJ/scripts/plot_maize_PB_groups_UpSet_02avn.py \
#    -t ${TAPPAS}/all_genotype_flag_tappas_DE_gene_results.csv \
#    -f ${IND}/flag_on_off_rsem_expression_isoforms_TPM_5_full_table.tsv \
#    -a ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_id_2_gene_id.csv \
#    -n ${REF}/B73_v4_geneIDs_geneNames_02amm.txt \
#    -v -m \
#    -p ${TAPPAS}/zmtr_fsm_ism_nic_nnc_consol

for TRT in Amb Ele; do
    python $PROJ/scripts/merge_GO_tappas_results.py \
        -t ${TAPPAS}/zmtr_fsm_ism_nic_nnc_consol_${TRT}_only.csv \
        --type DIU \
        --GOid ${REF}/sasdata/GO_IDs_catted.csv \
        --GOterm $PROJ/go_enrichment/GO_table.txt \
        --output-full ${TAPPAS}/zmtr_fsm_ism_nic_nnc_consol_${TRT}_only_w_GO.csv
done
