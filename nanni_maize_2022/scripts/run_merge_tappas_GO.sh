#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge tappAS output files with GO

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
TAPPAS=${IND}/tappas_output_TMM_norm_1CPMfilter
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info

#python $PROJ/scripts/merge_GO_tappas_results.py \
#    -t ${TAPPAS}/all_genotype_flag_tappas_DE_gene_results.csv \
#    --GOid ${REF}/RefGenV4/sasdata/GO_IDs_catted.csv \
#    --GOterm $PROJ/go_enrichment/GO_table.txt \
#    --output-full ${TAPPAS}/all_genotype_flag_tappas_DE_gene_results_w_GO.csv \
#    --output-5 ${TAPPAS}/zmtr_fsm_ism_nic_nnc_consol_DE_all5genotype_w_GO.csv

#python $PROJ/scripts/merge_GO_tappas_results.py \
#    -t ${TAPPAS}/all_genotype_flag_tappas_DIU_gene_results.csv \
#    --type DIU \
#    --GOid ${REF}/RefGenV4/sasdata/GO_IDs_catted.csv \
#    --GOterm $PROJ/go_enrichment/GO_table.txt \
#    --output-full ${TAPPAS}/all_genotype_flag_tappas_DIU_gene_results_w_GO.csv

python $PROJ/scripts/merge_GO_tappas_results.py \
    -t ${TAPPAS}/all_genotype_flag_tappas_DIU_gene_results_atLeast1geno_DIU_majorIsoformSwitch.csv \
    --type DIU \
    --GOid ${REF}/RefGenV4/sasdata/GO_IDs_catted.csv \
    --GOterm $PROJ/go_enrichment/GO_table.txt \
    --output-full ${TAPPAS}/all_genotype_flag_tappas_DIU_gene_results_atLeast1geno_DIU_majorIsoformSwitch_w_GO.csv
