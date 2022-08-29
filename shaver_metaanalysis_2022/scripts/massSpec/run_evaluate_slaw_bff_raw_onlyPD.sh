#!/bin/bash


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts


INPUT=$PROJ/SLAW_UGA_Output/bff_filtering
QC=$PROJ/SLAW_UGA_Output/qc_raw_${VAR}_onlyPD
    mkdir -p $QC

echo "

CV unlogged

"
python $SCRIPTS/coefficient_variation_flags.py \
    --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group newBatch \
    --CVcutoff 0.1 \
    --figure $QC/CV_figure_batch_${VAR}_PH_raw.pdf \
    --flag $QC/CV_flags_batch_${VAR}_PH_raw.tsv


    echo "

    CV unlogged

    "
    python $SCRIPTS/coefficient_variation_flags.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType \
        --CVcutoff 0.1 \
        --figure $QC/CV_figure_sampleType_${VAR}_PH_raw.pdf \
        --flag $QC/CV_flags_sampleType_${VAR}_PH_raw.tsv


    ## SED
    echo "

    SED batch

    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group batch \
        --figure $QC/SED_figure_batch_${VAR}_PH_raw.pdf \
        --SEDtoMean $QC/SEDtoMean_batch_${VAR}_PH_raw.tsv \
        --SEDpairwise $QC/SEDpairwise_batch_${VAR}_PH_raw.tsv \
        --per 0.8


    ## SED
    echo "

    SED sampleType
    
    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType \
        --figure $QC/SED_figure_sampleType_${VAR}_PH_raw.pdf \
        --SEDtoMean $QC/SEDtoMean_sampleType_${VAR}_PH_raw.tsv \
        --SEDpairwise $QC/SEDpairwise_sampleType_${VAR}_PH_raw.tsv \
        --per 0.8

    ## SED
    echo "

    SED sampleTypeBatch [wont run on pool, less than 3 per batch]
    
    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType_batch \
        --figure $QC/SED_figure_sampleTypeBatch_${VAR}_PH_raw.pdf \
        --SEDtoMean $QC/SEDtoMean_sampleTypeBatch_${VAR}_PH_raw.tsv \
        --SEDpairwise $QC/SEDpairwise_sampleTypeBatch_${VAR}_PH_raw.tsv \
        --per 0.8

    ## PCA
    echo "

    PCA sampleType
    
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType \
        --load_out $QC/PCA_loading_sampleType_${VAR}_PH_raw.tsv \
        --score_out $QC/PCA_scores_sampleType_${VAR}_PH_raw.tsv \
        --summary_out $QC/PCA_summary_sampleType_${VAR}_PH_raw.tsv \
        --figure $QC/PCA_figure_sampleType_${VAR}_PH_raw.pdf


  ## PCA
    echo "
    PCA batch
    
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group batch \
        --load_out $QC/PCA_loading_batch_${VAR}_PH_raw.tsv \
        --score_out $QC/PCA_scores_batch_${VAR}_PH_raw.tsv \
        --summary_out $QC/PCA_summary_batch_${VAR}_PH_raw.tsv \
        --figure $QC/PCA_figure_batch_${VAR}_PH_raw.pdf


  ## PCA
    echo "

    PCA sampleType_batch
    
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType_batch \
        --load_out $QC/PCA_loading_sampleTypeBatch_${VAR}_PH_raw.tsv \
        --score_out $QC/PCA_scores_sampleTypeBatch_${VAR}_PH_raw.tsv \
        --summary_out $QC/PCA_summary_sampleTypeBatch_${VAR}_PH_raw.tsv \
        --figure $QC/PCA_figure_sampleTypeBatch_${VAR}_PH_raw.pdf

  ## distribution samples 
    echo "

    Distribution samples rank sampleType
    
    "
    python $SCRIPTS/distribution_samples.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType \
        --figure $QC/distSamples_sampleType_${VAR}_PH_raw.pdf

  ## distribution features
    echo "
 
    Distribution features rank sampleType
    
    "
    python $SCRIPTS/distribution_features.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group sampleType \
        --figure $QC/distFeatures_sampleType_${VAR}_PH_raw.pdf


  ## distribution samples 
    echo "

    Distribution samples batch
    
    "
    python $SCRIPTS/distribution_samples.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group newBatch \
        --figure $QC/distSamples_batch_${VAR}_PH_raw.pdf

  ## distribution features
    echo "
 
    Distribution features batch
    
    "
    python $SCRIPTS/distribution_features.py \
        --input $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
        --design $DESIGN \
        --ID UniqueID \
        --group newBatch \
        --figure $QC/distFeatures_batch_${VAR}_PH_raw.pdf

