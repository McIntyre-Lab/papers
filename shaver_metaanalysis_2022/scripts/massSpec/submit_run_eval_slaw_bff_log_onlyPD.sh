#!/bin/bash


##QC on only PD1074 and PD1074 pools

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
LOG=$PROJ/scripts/logs_slaw
    mkdir -p $LOG


## RP NEG:
export VARDIR=RP_NEG
export VAR=rp_neg
export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_NEG_slaw.tsv

## run QC raw 
source ${PROJ}/scripts/massSpec/run_evaluate_slaw_bff_log_onlyPD.sh 2>&1 > $LOG/rp_neg_log_qc.log

## RP POS
export VARDIR=RP_POS
export VAR=rp_pos
export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw.tsv

## run QC raw noBlanks
source ${PROJ}/scripts/massSpec/run_evaluate_slaw_bff_log_onlyPD.sh 2>&1 > $LOG/rp_pos_log_qc.log

## HILIC POS
export VARDIR=HILIC_POS
export VAR=hilic_pos
export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/hilic_gtech/dsgn_GT_HILIC_POS_slaw.tsv

## run QC raw noBlanks
source ${PROJ}/scripts/massSpec/run_evaluate_slaw_bff_log_onlyPD.sh 2>&1 > $LOG/hilic_pos_log_qc.log

## HILIC NEG
#export VARDIR=HILIC_NEG
#export VAR=hilic_neg
#export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/hilic_gtech/dsgn_GT_HILIC_NEG_slaw.tsv

## run QC raw noBlanks
#source ${PROJ}/scripts/massSpec/run_evaluate_slaw_bff_log_onlyPD.sh 2>&1 > $LOG/hilic_neg_log_qc.log

