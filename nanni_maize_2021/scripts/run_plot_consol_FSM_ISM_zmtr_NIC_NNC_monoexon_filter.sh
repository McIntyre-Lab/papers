#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Plot transcripts per gene, number of exon regions per gene,
##     number of exonic fragments per gene, and proportion of
##     varying exonic regions/fragments per gene for consolidated
##     representative reference transcripts of FSM/ISM
##     and the monoexon filtered NIC/NNC

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=$PROJ/scripts
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
IND=$PROJ/sqanti_classification_category_subset/EA_annotations
OUTD=$PROJ/sqanti_classification_category_subset/plot_consol_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
    mkdir -p ${OUTD}

## Subset EA annotations for consolidated transcriptome
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/zmtr_fsm_ism_nic_nnc_event2transcript2gene_index.csv \
#    -v event_id \
#    -l ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_list_not_retained.csv \
#    -o ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/zmtr_fsm_ism_nic_nnc_fusion_annotations.csv \
#    -v fusion_id \
#    -l ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_list_not_retained.csv \
#    -o ${IND}/zmtr_fsm_ism_nic_nnc_consol_fusion_annotations.csv
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/zmtr_fsm_ism_nic_nnc_exon_fragment_annotations.csv \
#    -v fragment_id \
#    -l ${DIST}/zmtr_fsm_ism_nic_nnc_transcript_list_not_retained.csv \
#    -o ${IND}/zmtr_fsm_ism_nic_nnc_consol_exon_fragment_annotations.csv

## Print total transcripts and genes represented
##     NIC/NNC transcripts have "PB" in the name,
##     FSM/ISM representative reference transcripts do not
echo "
$(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv | \
    sed s/'|'/'\n'/g | sort | uniq | wc -l) transcripts represented ($(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv | \
    sed s/'|'/'\n'/g | sort | uniq | grep -v "PB" | wc -l) consolidated FSM/ISM representative reference transcripts and $(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv | \
    sed s/'|'/'\n'/g | sort | uniq | grep "PB" | wc -l) consolidated filtered NIC/NNC)
$(awk -F "," 'NR!=1{print $3}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv | \
    sed s/'|'/'\n'/g | sort | uniq | wc -l) reference genes represented
"

## Plot
python ${SCRIPTS}/EA_histgram_visualize_07lzh.py \
    -ig ${IND}/zmtr_fsm_ism_nic_nnc_consol_event2transcript2gene_index.csv \
    -ie ${IND}/zmtr_fsm_ism_nic_nnc_consol_fusion_annotations.csv \
    -if ${IND}/zmtr_fsm_ism_nic_nnc_consol_exon_fragment_annotations.csv \
    --prefix zmtr_fsm_ism_nic_nnc_consol \
    -o ${OUTD} \
    -s -l -r
