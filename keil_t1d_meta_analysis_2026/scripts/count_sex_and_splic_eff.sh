#!/bin/bash

##Extract pvalues from meta-analysis summary
PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

##Set paths
IND=${PROJ}/quantify_t1d_pacbio_transcripts
SCRIPTS=${PROJ}/scripts


##cp files for extraction to roz
for CELL in CD4 CD8; do
    python ${SCRIPTS}/count_sex_and_splicing_effect.py \
        -i ${IND}/metaanalysis_summary_${CELL}.csv \
        -o ${IND} \
        -p ${CELL}

done
