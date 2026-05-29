#!/bin/bash

###  run DE model
    ## fragments!!

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

#PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType
SASP=${PROJ}/sas_programs

#PROGRAM=model_frag_DE_age_sex_04amm
#sas -log ${PROJ}/sas_logs/${PROGRAM}.log \
#    -work /sas/SASWORK \
#    -sysin ${SASP}/${PROGRAM}.sas


PROGRAM2=prep_pdiffs_4_meta_analysis_05amm
sas -log ${PROJ}/sas_logs/${PROGRAM2}.log \
    -work /sas/SASWORK \
    -sysin ${SASP}/${PROGRAM2}.sas

