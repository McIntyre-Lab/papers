#!/bin/bash


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts


INPUT=$PROJ/SLAW_UGA_Output/bff_filtering
QC=$PROJ/SLAW_UGA_Output/qc_log_${VAR}
    mkdir -p $QC


echo "

log transform

"
python $SCRIPTS/log_and_glog_transformation.py \
    --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
    --design $DESIGN/dsgn_GT_${VARDIR}_slaw.tsv \
    --uniqID UniqueID \
    --transformation log \
    --log_base log2 \
    --oname $INPUT/BFF_corrected_${VAR}_PH_log.tsv


for i in batch strain_w_poolPD
do

if [[ $i = batch ]]; then 
    DSGN=dsgn_GT_${VARDIR}_slaw.tsv
else 
    DSGN=dsgn_GT_${VARDIR}_slaw_strainWPD.tsv
fi


echo "

CV by ${i}

"
python $SCRIPTS/coefficient_variation_flags.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_log.tsv \
    --design $DESIGN/$DSGN \
    --ID UniqueID \
    --group ${i} \
    --CVcutoff 0.1 \
    --figure $QC/CV_figure_${i}_${VAR}_PH_log.pdf \
    --flag $QC/CV_flags_${i}_${VAR}_PH_log.tsv



    ## SED
    echo "

    SED by ${i}

    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/BFF_corrected_${VAR}_PH_log.tsv \
        --design $DESIGN/$DSGN \
        --ID UniqueID \
        --group ${i} \
        --figure $QC/SED_figure_${i}_${VAR}_PH_log.pdf \
        --SEDtoMean $QC/SEDtoMean_${i}_${VAR}_PH_log.tsv \
        --SEDpairwise $QC/SEDpairwise_batch_${VAR}_PH_log.tsv \
        --per 0.8

## PCA
    echo "

    PCA by ${i}
    
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/BFF_corrected_${VAR}_PH_log.tsv \
        --design $DESIGN/$DSGN \
        --ID UniqueID \
        --group ${i} \
        --load_out $QC/PCA_loading_${i}_${VAR}_PH_log.tsv \
        --score_out $QC/PCA_scores_${i}_${VAR}_PH_log.tsv \
        --summary_out $QC/PCA_summary_${i}_${VAR}_PH_log.tsv \
        --figure $QC/PCA_figure_${i}_${VAR}_PH_log.pdf



  ## distribution samples 
    echo "

    Distribution samples  by ${i}
    
    "
    python $SCRIPTS/distribution_samples.py \
        --input $INPUT/BFF_corrected_${VAR}_PH_log.tsv \
        --design $DESIGN/$DSGN \
        --ID UniqueID \
        --group sampleType \
        --figure $QC/distSamples_${i}_${VAR}_PH_log.pdf

  ## distribution features
    echo "
 
    Distribution features by ${i}
    
    "
    python $SCRIPTS/distribution_features.py \
        --input $INPUT/BFF_corrected_${VAR}_PH_log.tsv \
        --design $DESIGN/$DSGN \
        --ID UniqueID \
        --group sampleType \
        --figure $QC/distFeatures_${i}_${VAR}_PH_log.pdf

done