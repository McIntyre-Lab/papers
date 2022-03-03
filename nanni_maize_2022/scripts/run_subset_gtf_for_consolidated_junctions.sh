#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Subset GTF file for consolidated transcriptome
## *** Will be used for RSEM quantification

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=$PROJ/scripts
IND=$PROJ/sqanti_classification_category_subset
DIST=${IND}/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter

## Subset GTF for consolidated transcriptome
#python ${SCRIPTS}/subset_gtf_simple_02avn.py \
#    -g ${IND}/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.gtf \
#    -t transcript_id \
#    -e ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_list_not_retained.csv \
#    -o ${IND}/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf

#rsync -ulvz ${IND}/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf \
#    adalena.nanni@hpg.rc.ufl.edu:/blue/mcintyre/share/maize_ainsworth/.

## Count transcripts/gene for consolidated transcriptome
python ${SCRIPTS}/count_GTF_transcript_per_gene.py \
    -g ${IND}/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf \
    -o ${IND}/consol_transcriptome_transcript_per_gene_freq.txt
