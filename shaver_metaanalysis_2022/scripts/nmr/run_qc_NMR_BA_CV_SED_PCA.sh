#!/bin/bash


#bash
#source activate /TB14/TB14/gait_gm/galaxy/database/dependencies/_conda/envs/__gaitGm@1.0.1
export PATH=/TB14/TB14/conda_envs/gait-gm/bin:$PATH

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16
SCRIPTS=/TB14/TB14/github/SECIMTools/src/scripts


#BA, SED, CV, PCA

INPUT=$PROJ/nmr
DESIGN=$INPUT/design_files_nmr
OUTPUT=$INPUT/qc
    mkdir -p $OUTPUT

for i in CDCL3 D2O
do

: <<'END'
echo "log transform ${i}
"
python $SCRIPTS/log_and_glog_transformation.py \
    --input $INPUT/ph_NMR_${i}_shortNames.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --uniqID ppm \
    --transformation log \
    --log_base log2 \
    --oname $INPUT/ph_NMR_${i}_shortNames_log.tsv

## run bland altman plot
echo "BA by genotype $i
"
python $SCRIPTS/bland_altman_plot.py \
    --input $INPUT/ph_NMR_${i}_shortNames_log.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --figure $OUTPUT/BA_plots_genotype_ph_NMR_${i}_log.pdf \
    --flag_dist $OUTPUT/BA_dist_flags_genotype_ph_NMR_${i}_log.tsv \
    --flag_sample $OUTPUT/BA_sample_flags_genotype_ph_NMR_${i}_log.tsv \
    --flag_feature $OUTPUT/BA_feature_flags_genotype_ph_NMR_${i}_log.tsv \
    --prop_feature $OUTPUT/BA_propFeature_genotype_ph_NMR_${i}_log.tsv \
    --prop_sample $OUTPUT/BA_propSample_genotype_ph_NMR_${i}_log.tsv

## CV
echo "CV ${i} group = genotype
"
python $SCRIPTS/coefficient_variation_flags.py \
    --input $INPUT/ph_NMR_${i}_shortNames_log.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --CVcutoff 0.1 \
    --figure $OUTPUT/CV_figure_ph_NMR_${i}.pdf \
    --flag $OUTPUT/CV_flags_ph_NMR_${i}.tsv 

## SED
    echo "SED ${i} order = run_order
    "
    python $SCRIPTS/standardized_euclidean_distance.py \
        --input $INPUT/ph_NMR_${i}_shortNames.tsv  \
        --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
        --ID ppm \
        --group genotype \
        --order run_order \
        --figure $OUTPUT/SED_figure_ph_NMR_${i}.pdf \
        --SEDtoMean $OUTPUT/SEDtoMean_ph_NMR_${i}.tsv \
        --SEDpairwise $OUTPUT/SEDpairwise_ph_NMR_${i}.tsv \
        --per 0.8
END

## PCA
    echo "PCA $i
    "
    python $SCRIPTS/principal_component_analysis.py \
        --input $INPUT/ph_NMR_${i}_shortNames.tsv \
        --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
        --ID ppm \
        --group genotype \
        --load_out $OUTPUT/PCA_loading_ph_NMR_${i}.tsv \
        --score_out $OUTPUT/PCA_scores_ph_NMR_${i}.tsv \
        --summary_out $OUTPUT/PCA_summary_ph_NMR_${i}.tsv \
        --figure $OUTPUT/PCA_figure_pos_ph_NMR_${i}.pdf

: <<'END'

## distribution samples 
echo "Distribution samples $i
"
python $SCRIPTS/distribution_samples.py \
    --input $INPUT/ph_NMR_${i}_shortNames.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --order run_order \
    --figure $OUTPUT/distSamples_ph_NMR_${i}.pdf

## distribution features
echo "Distribution features $i
"
python $SCRIPTS/distribution_features.py \
    --input $INPUT/ph_NMR_${i}_shortNames.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --figure $OUTPUT/distFeatures_ph_NMR_${i}.pdf

echo "Distribution samples log $i
"
python $SCRIPTS/distribution_samples.py \
    --input $INPUT/ph_NMR_${i}_shortNames_log.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --order run_order \
    --figure $OUTPUT/distSamples_ph_NMR_${i}_log.pdf

## distribution features
echo "Distribution features log $i
"
python $SCRIPTS/distribution_features.py \
    --input $INPUT/ph_NMR_${i}_shortNames_log.tsv \
    --design $DESIGN/dsgn_NMR_${i}_shortNames.txt \
    --ID ppm \
    --group genotype \
    --figure $OUTPUT/distFeatures_ph_NMR_${i}_log.pdf

END

done
