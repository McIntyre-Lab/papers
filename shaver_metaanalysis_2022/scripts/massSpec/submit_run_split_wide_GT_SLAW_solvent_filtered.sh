#!/bin/bash


##
PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16

## RP NEG:
export VAR=rp_neg

## split starting datasets
source ${PROJ}/scripts/massSpec/run_split_wide_GT_SLAW_solvent_filtered.sh

## RP POS
export VAR=rp_pos

## split starting datasets
source ${PROJ}/scripts/massSpec/run_split_wide_GT_SLAW_solvent_filtered.sh

## HILIC POS
export VAR=hilic_pos

## split starting datasets
source ${PROJ}/scripts/massSpec/run_split_wide_GT_SLAW_solvent_filtered.sh

## HILIC NEG
export VAR=hilic_neg

## split starting datasets
source ${PROJ}/scripts/massSpec/run_split_wide_GT_SLAW_solvent_filtered.sh

