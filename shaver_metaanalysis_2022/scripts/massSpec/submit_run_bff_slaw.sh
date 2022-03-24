#!/bin/bash



PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
LOG=$PROJ/scripts/logs_slaw
    mkdir -p $LOG


## RP NEG:
#export VAR=rp_neg
#export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_NEG_slaw.tsv

#source ${PROJ}/scripts/massSpec/run_bff_slaw.sh 2>&1 > $LOG/rp_neg_bff.log


## RP POS
#export VAR=rp_pos
#export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/dsgn_GT_RP_POS_slaw.tsv

#source ${PROJ}/scripts/massSpec/run_bff_slaw.sh 2>&1 > $LOG/rp_pos_bff.log


## HILIC POS
#export VAR=hilic_pos
#export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/hilic_gtech/dsgn_GT_HILIC_POS_slaw.tsv

#source ${PROJ}/scripts/massSpec/run_bff_slaw.sh 2>&1 > $LOG/hilic_pos_bff.log

## HILIC NEG
export VAR=hilic_neg
export DESIGN=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/design_files/hilic_gtech/dsgn_GT_HILIC_NEG_slaw.tsv

## run QC ranked noBlanks
source ${PROJ}/scripts/massSpec/run_bff_slaw.sh 2>&1 > $LOG/hilic_neg_bff.log

