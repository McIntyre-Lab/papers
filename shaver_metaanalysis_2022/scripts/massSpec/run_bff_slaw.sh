#!/bin/bash

## design and var specified in submission script

#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

INPUT=$PROJ/SLAW_UGA_Output/text_data
QC=$PROJ/SLAW_UGA_Output/bff_filtering
    mkdir -p $QC
ROZ=$PROJ/SLAW_UGA_Output/roz_amm
    mkdir -p $ROZ

echo "

dump library std, solera, batch_pool from all starting datasets
drop aos_41 and aos53 from RP neg and pos datasets

"
if [[ $VAR = *rp* ]] ; then
   DROPS=`egrep -i 'aos_41|aos_53|solera|batch|std' $DESIGN | cut -f1`
else
   DROPS=`egrep -i 'solera|batch|std' $DESIGN | cut -f1`
fi

python $SCRIPTS/remove_user_specified_row_col.py \
    --input $INPUT/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_sbys.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --row \
    --col $DROPS \
    --outWide $ROZ/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_sbys_data.tsv

echo "

Run BFF

    "
## BFF with default bff value of 5000 and default criterion value of 100
    ## outs is output for sample only (blank samples removed) file
python $SCRIPTS/blank_feature_filtering_flags_with_LOD.py \
    --input $ROZ/meta_analysis_${VAR}_SLAW_output_0xBFF_PH_sbys_data.tsv \
    --design $DESIGN \
    --uniqID UniqueID \
    --group sampleType \
    --blank blank \
    --bff 5000 \
    --criteria 100 \
    --outflags $QC/BFF_flag_file_${VAR}.tsv \
    --outbff $QC/BFF_file_${VAR}.tsv \
    --outlod $QC/BFF_lod_file_${VAR}.tsv \
    --outs $QC/${VAR}_PH_sbys_noBlank.tsv


## drop rows where flag bff_PD1074_off = 1 or flag_bff_pool_pd1074_off = 1 
  ## condition of drop :
        ## condition = 0 equal to
        ## condition = 1 greater than
        ## condition = 2 less than
  ## cutoff value - rows with flag_value equal to (0), greater than (1) or less than (2) the cutoff value will be dropped

echo "

drop selected features -- drop if flag_bff_PD1074_off = 1 or flag_bff_pool_pd1074_off = 1

    "
python $SCRIPTS/remove_selected_features_samples.py \
    --input $QC/${VAR}_PH_sbys_noBlank.tsv \
    --design $DESIGN \
    --flags $QC/BFF_flag_file_${VAR}.tsv \
    --ID UniqueID \
    --flagUniqID UniqueID  \
    --flagDrop flag_bff_PD1074_off \
    --flagfiletype row \
    --value 1 \
    --condition 0 \
    --outWide $QC/BFF_corrected_${VAR}_PD.tsv \
    --outFlags $QC/flags_BFF_corrected_${VAR}_PD.tsv


python $SCRIPTS/remove_selected_features_samples.py \
    --input $QC/BFF_corrected_${VAR}_PD.tsv \
    --design $DESIGN \
    --flags $QC/BFF_flag_file_${VAR}.tsv \
    --ID UniqueID \
    --flagUniqID UniqueID  \
    --flagDrop flag_bff_pool_pd1074_off \
    --flagfiletype row \
    --value 1 \
    --condition 0 \
    --outWide $QC/BFF_corrected_${VAR}_PD_poolPD.tsv \
    --outFlags $QC/flags_BFF_corrected_${VAR}_PD_poolPD.tsv

rm -r $ROZ
