#!/bin/bash


#export Path to Secim conda env 
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__gait-gm@21.7.22/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/galaxy_27oct21/galaxy/tools/SECIMTools/src/scripts
GAIT=/TB14/TB14/galaxy_27oct21/galaxy/tools/gait-gm/src/scripts

INPUT=$PROJ/SLAW_UGA_Output/bff_filtering
QC=$PROJ/SLAW_UGA_Output/qc_${VAR}_rank_onlyPD
    mkdir -p $QC

echo "

rank transform bin 100

"
python $SCRIPTS/add_group_rank.py \
    --wide $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
    --design $DESIGN \
    --uniqID UniqueID \
    --ngroup 100 \
    --out $INPUT/BFF_corrected_${VAR}_PH_rank100.tsv


echo "

rank transform

"
python $SCRIPTS/add_group_rank.py \
    --wide $INPUT/BFF_corrected_${VAR}_PD_poolPD.tsv \
    --design $DESIGN \
    --uniqID UniqueID \
    --out $INPUT/BFF_corrected_${VAR}_PH_rank.tsv

echo "

CV

"
python $SCRIPTS/coefficient_variation_flags.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --CVcutoff 0.1 \
    --figure $QC/CV_figure_sampleType_${VAR}_PH_rank.pdf \
    --flag $QC/CV_flags_sampleType_${VAR}_PH_rank.tsv

echo "

CV

"
python $SCRIPTS/coefficient_variation_flags.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group newBatch \
    --CVcutoff 0.1 \
    --figure $QC/CV_figure_batch_${VAR}_PH_rank.pdf \
    --flag $QC/CV_flags_batch_${VAR}_PH_rank.tsv

## SED
    echo "

SED batch

"
python $SCRIPTS/standardized_euclidean_distance.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group batch \
    --figure $QC/SED_figure_batch_${VAR}_PH_rank.pdf \
    --SEDtoMean $QC/SEDtoMean_batch_${VAR}_PH_rank.tsv \
    --SEDpairwise $QC/SEDpairwise_batch_${VAR}_PH_rank.tsv \
    --per 0.8

## SED
echo "

SED sampleType

"
python $SCRIPTS/standardized_euclidean_distance.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --figure $QC/SED_figure_sampleType_${VAR}_PH_rank.pdf \
    --SEDtoMean $QC/SEDtoMean_sampleType_${VAR}_PH_rank.tsv \
    --SEDpairwise $QC/SEDpairwise_sampleType_${VAR}_PH_rank.tsv \
    --per 0.8

## SED
echo "

SED sampleType_batch

"
python $SCRIPTS/standardized_euclidean_distance.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType_batch \
    --figure $QC/SED_figure_sampleTypeBatch_${VAR}_PH_rank.pdf \
    --SEDtoMean $QC/SEDtoMean_sampleTypeBatch_${VAR}_PH_rank.tsv \
    --SEDpairwise $QC/SEDpairwise_sampleTypeBatch_${VAR}_PH_rank.tsv \
    --per 0.8

## PCA
echo "

PCA sampleType

"
python $SCRIPTS/principal_component_analysis.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --load_out $QC/PCA_loading_sampleType_${VAR}_PH_rank.tsv \
    --score_out $QC/PCA_scores_sampleType_${VAR}_PH_rank.tsv \
    --summary_out $QC/PCA_summary_sampleType_${VAR}_PH_rank.tsv \
    --figure $QC/PCA_figure_sampleType_${VAR}_PH_rank.pdf

## PCA
echo "

PCA batch

"
 python $SCRIPTS/principal_component_analysis.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group batch \
    --load_out $QC/PCA_loading_batch_${VAR}_PH_rank.tsv \
    --score_out $QC/PCA_scores_batch_${VAR}_PH_rank.tsv \
    --summary_out $QC/PCA_summary_batch_${VAR}_PH_rank.tsv \
    --figure $QC/PCA_figure_batch_${VAR}_PH_rank.pdf

## PCA
echo "

PCA sampleType_batch

"
python $SCRIPTS/principal_component_analysis.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType_batch \
    --load_out $QC/PCA_loading_sampleTypeBatch_${VAR}_PH_rank.tsv \
    --score_out $QC/PCA_scores_sampleTypeBatch_${VAR}_PH_rank.tsv \
    --summary_out $QC/PCA_summary_sampleTypeBatch_${VAR}_PH_rank.tsv \
    --figure $QC/PCA_figure_sampleTypeBatch_${VAR}_PH_rank.pdf


 ## distribution samples 
 echo "

Distribution samples rank sampleType

"
python $SCRIPTS/distribution_samples.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --figure $QC/distSamples_sampleType_${VAR}_PH_rank.pdf

## distribution features
echo "

Distribution features rank sampleType

"
python $SCRIPTS/distribution_features.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --figure $QC/distFeatures_sampleType_${VAR}_PH_rank.pdf

## distribution samples 
echo "

Distribution samples rank sampleType

"
python $SCRIPTS/distribution_samples.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group newBatch \
    --figure $QC/distSamples_batch_${VAR}_PH_rank.pdf

## distribution features
echo "

Distribution features rank sampleType

"
python $SCRIPTS/distribution_features.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group newBatch \
    --figure $QC/distFeatures_batch_${VAR}_PH_rank.pdf

## switch to secim env for BA plots of PD1074
export PATH=/TB14/TB14/galaxy_27oct21/galaxy/database/dependencies/_conda/envs/__secimtools@21.2.21:$PATH

echo "

BA plots

"
python $SCRIPTS/bland_altman_plot.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank100.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --process_only pool_pd1074 \
    --figure $QC/BA_fig_sampleType_pd1074Pool_${VAR}_PH_rank100.pdf \
    --flag_dist $QC/BA_dist_flags_sampleType_pd1074Pool_${VAR}_PH_rank100.tsv \
    --flag_sample $QC/BA_sample_flags_sampleType_pd1074Pool_${VAR}_PH_rank100.tsv \
    --flag_feature $QC/BA_feature_flags_sampleType_pd1074Pool_${VAR}_PH_rank100.tsv \
    --prop_feature $QC/BA_propFeature_sampleType_pd1074Pool_${VAR}_PH_rank100.tsv \
    --prop_sample $QC/BA_propSample_sampleType_pd1074Pool_${VAR}_PH_rank100.tsv

 python $SCRIPTS/bland_altman_plot.py \
    --input $INPUT/BFF_corrected_${VAR}_PH_rank100.tsv \
    --design $DESIGN \
    --ID UniqueID \
    --group sampleType \
    --process_only PD1074 \
    --figure $QC/BA_fig_sampleType_pd1074_${VAR}_PH_rank100.pdf \
    --flag_dist $QC/BA_dist_flags_sampleType_pd1074_${VAR}_PH_rank100.tsv \
    --flag_sample $QC/BA_sample_flags_sampleType_pd1074_${VAR}_PH_rank100.tsv \
    --flag_feature $QC/BA_feature_flags_sampleType_pd1074_${VAR}_PH_rank100.tsv \
    --prop_feature $QC/BA_propFeature_sampleType_pd1074_${VAR}_PH_rank100.tsv \
    --prop_sample $QC/BA_propSample_sampleType_pd1074_${VAR}_PH_rank100.tsv

