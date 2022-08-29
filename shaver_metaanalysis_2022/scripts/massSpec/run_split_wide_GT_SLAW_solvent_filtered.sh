#!/bin/bash


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH


PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

DESIGN=$PROJ/design_files

INPUT=$PROJ/SLAW_UGA_Output/filtered_datasets
OUTPUT=$PROJ/SLAW_UGA_Output/filtered_datasets_processed
    mkdir -p $OUTPUT

ROZ=$PROJ/SLAW_UGA_Output/roz_amm
    mkdir -p $ROZ

### csv to tsv
sed 's@,@\t@g' $INPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt.txt > $INPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt.tsv


## figure out number of data columns
FIRST=10
LAST=`awk '{print NF; exit}' $INPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt.tsv`

seq $FIRST $LAST > $ROZ/range_${VAR}.tsv | sed -z 's@\n@,@g;s@,$@\n@' $ROZ/range_${VAR}.tsv > $ROZ/range_${FILE}.csv
RANGE=`cat $ROZ/range_${FILE}.csv`

## create wide dataset, anno (non-data columns) and design file start
python ${GAIT}/split_wide_dataset.py \
    --input ${INPUT}/BFF_corrected_${VAR}_with_rt_mz_solv_filt.tsv \
    --samples $RANGE \
    --ID uniqueID \
    --wide $OUTPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt_sbys.tsv \
    --design $OUTPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt_design.tsv \
    --annot $OUTPUT/BFF_corrected_${VAR}_with_rt_mz_solv_filt_annot.tsv

rm -r $ROZ
