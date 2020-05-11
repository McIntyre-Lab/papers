#!/bin/bash



## summarize snp counts for matching and non-matching isolate pairs

PROJ=/home/ammorse/mclab/SHARE/McIntyre_Lab/staph/seoung_ho_project
SASP=$PROJ/sas_programs
LOGS=$PROJ/sas_logs
    mkdir -p $LOGS


PRG=summarize_snp_counts_07amm
sas -log ${LOGS}/${PRG}.log -work /home/ammorse/SASWORK -sysin ${SASP}/${PRG}.sas
