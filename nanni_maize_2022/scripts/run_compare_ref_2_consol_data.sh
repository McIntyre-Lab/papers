#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Plot transcripts per gene comparisons

### Set Directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/sqanti_classification_category_subset/plot_consol_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
INEA=$PROJ/sqanti_classification_category_subset/EA_annotations
WANGEA=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/maize_leaf_150bp_annotations

## Get titles
DATANAME="Maize Leaf (5 genotypes, 2 treatments)"
REFNAME="Maize B73 Reference"
WANGNAME="Wang, et al. 2018 Maize B73 Leaf"

## Get event files
REFEVENT=~/mclab/SHARE/McIntyre_Lab/useful_maize_info/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_event2transcript2gene_index.csv
DATAEVENT=${INEA}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv
WANGEVENT=${WANGEA}/maize_leaf_consol_event2transcript2gene_index.csv

## Plot transcripts per gene for consolidated transcriptome vs. reference
echo "
${REFNAME} vs. ${DATANAME}"
python ${SCRIPTS}/plot_compare_xcrpt_per_gene_02avn.py \
    -1 ${REFEVENT} \
    -2 ${DATAEVENT} \
    -n1 "${REFNAME}" \
    -n2 "${DATANAME}" \
    --hist-1-2 \
    --scatter \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_consol_vs_reference

## Plot transcripts per gene for consolidated transcriptome vs. 
##     Wang, et al 2018 maize B73 leaf transcriptome (representative reference
##     transcripts of FSM/ISM)
echo "
${WANGNAME} vs. ${DATANAME}"
python ${SCRIPTS}/plot_compare_xcrpt_per_gene_02avn.py \
    -1 ${WANGEVENT} \
    -2 ${DATAEVENT} \
    -n1 "${WANGNAME}" \
    -n2 "${DATANAME}" \
    --hist-1-2 \
    --hist-2-1 \
    --scatter \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_consol_vs_MaizeWang2018_leaf
