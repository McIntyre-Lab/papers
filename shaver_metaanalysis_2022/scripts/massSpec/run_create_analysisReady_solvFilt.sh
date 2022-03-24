#!/bin/bash

## create datasets ready for meta analysis
##    RP neg, RP pos, HILIC pos
##    NOT including HILIC neg for now

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

INPUT=$PROJ/SLAW_UGA_Output/filtered_datasets_processed
READY=$PROJ/SLAW_UGA_Output/analysisReady_datasets
    mkdir -p $READY


if [[ $VAR = *hilic_pos* ]] ; then
   DROPS=`egrep -i 'b1_aos_25|b5_aos_49' $DESIGN/dsgn_GT_${VARDIR}_slaw.tsv | cut -f1`
else
   DROPS=`egrep -i 'b5_aos_49' $DESIGN/dsgn_GT_${VARDIR}_slaw.tsv | cut -f1`
fi

python $SCRIPTS/remove_user_specified_row_col.py \
    --input $INPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt_sbys.tsv \
    --design $DESIGN/dsgn_GT_${VARDIR}_slaw.tsv \
    --ID uniqueID \
    --row \
    --col $DROPS \
    --outWide $READY/${VAR}_analysisReady_sbys.tsv


python $SCRIPTS/add_group_rank.py \
    --wide $READY/${VAR}_analysisReady_sbys.tsv \
    --design $DESIGN/dsgn_GT_${VARDIR}_slaw.tsv \
    --uniqID uniqueID \
    --out $READY/${VAR}_analysisReady_rank_sbys.tsv
