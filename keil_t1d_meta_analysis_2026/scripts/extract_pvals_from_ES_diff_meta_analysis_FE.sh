#!/bin/bash

##Extract pvalues from meta-analysis summary
PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

##Set paths - overall metanalysis
IND=${PROJ}/quantify_t1d_pacbio_transcripts/meta_analysis_sex_diff_summary
SCRIPTS=${PROJ}/scripts

ROZ=${IND}/roz
   mkdir -p ${ROZ}

##cp files for extraction to roz
for file in ${IND}/summary_ES_diff*; do
    cp "$file" ${ROZ}
done

#cp ${IND_ALL}/summary_sexMod_* ${ROZ_ALL}


#Extract pvalues and output CSV - overall meta
python ${SCRIPTS}/extract_pvals_from_meta_analysis_summary_no_moderator_02ndk.py \
    -d ${ROZ}\
    -s FE.txt\
    -x 0.05 \
    -o ${IND}/meta_analysis_pvals_ES_diff_heterogenity.csv

rm -r ${ROZ}


