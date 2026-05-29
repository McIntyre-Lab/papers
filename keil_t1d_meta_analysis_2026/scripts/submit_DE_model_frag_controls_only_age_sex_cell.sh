#!/bin/bash

###  run DE model
    ## fragments!!
    ## controls only, CD4 + CD8

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

#PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType

SASP=${PROJ}/sas_programs


PROGRAM=model_frag_DE_controls_only_age_sex_cell

sas -log ${PROJ}/sas_logs/${PROGRAM}.log \
    -work /sas/SASWORK \
    -sysin ${SASP}/${PROGRAM}.sas

