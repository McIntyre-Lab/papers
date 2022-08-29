#!/bin/bash


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH


PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

DESIGN=$PROJ/design_files

INPUT=$PROJ/SLAW_UGA_Output/original_data
OUTPUT=$PROJ/SLAW_UGA_Output/text_data
    mkdir -p $OUTPUT

ROZ=$PROJ/SLAW_UGA_Output/roz_amm
    mkdir -p $ROZ

### txt to tsv
cp $INPUT/meta_analysis_${VAR}_SLAW_output_0xBFF.txt $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF.tsv

## convert all lowercase
tr '[:upper:]' '[:lower:]' < $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF.tsv > $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_low.tsv

## figure out number of data columns
FIRST=37
LAST=`awk '{print NF; exit}' $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_low.tsv`

seq $FIRST $LAST > $ROZ/range_${VAR}.tsv | sed -z 's@\n@,@g;s@,$@\n@' $ROZ/range_${VAR}.tsv > $ROZ/range_${FILE}.csv
RANGE=`cat $ROZ/range_${FILE}.csv`

## create wide dataset, anno (non-data columns) and design file start
python ${GAIT}/split_wide_dataset.py \
    --input ${OUTPUT}/meta_analysis_${VAR}_SLAW_output_0xBFF_low.tsv \
    --samples $RANGE \
    --prefix met \
    --wide $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_sbys.tsv \
    --design $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_design.tsv \
    --annot $OUTPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_annot.tsv

rm -r $ROZ
