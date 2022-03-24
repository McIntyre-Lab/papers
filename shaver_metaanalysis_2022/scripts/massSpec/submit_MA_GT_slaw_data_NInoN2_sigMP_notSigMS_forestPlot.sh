#!/bin/bash


##  submission script to run MA on RP pos, RP neg, Hilic pos and Hilic neg

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
LOG=$PROJ/SLAW_UGA_Output/ma_logs
    mkdir -p $LOG

## RP POS:
export VARDIR=RP_POS
export VAR=rp_pos

source ${PROJ}/scripts/massSpec/meta_analysis/run_MA_slaw_by_pathway_FE_rank_massSpec_NInoN2_sigMP_notSigMS_forestPlot.sh  2>&1 > $LOG/rp_pos_ma_path_NInoN2_forest.log

## RP NEG:
export VARDIR=RP_NEG
export VAR=rp_neg


source ${PROJ}/scripts/massSpec/meta_analysis/run_MA_slaw_by_pathway_FE_rank_massSpec_NInoN2_sigMP_notSigMS_forestPlot.sh 2>&1 > $LOG/rp_pos_ma_path_NInoN2_forest.log

## HILIC POS:
export VARDIR=HILIC_POS
export VAR=hilic_pos

source ${PROJ}/scripts/massSpec/meta_analysis/run_MA_slaw_by_pathway_FE_rank_massSpec_NInoN2_sigMP_notSigMS_forestPlot.sh  2>&1 > $LOG/rp_pos_ma_path_NInoN2_forest.log
