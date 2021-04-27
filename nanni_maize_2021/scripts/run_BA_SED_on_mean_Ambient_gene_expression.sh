#!/bin/bash

## run BA and SED on meand ambient values

PROJ=/home/ammorse/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=/home/ammorse/TB14/galaxy/tools/GalaxyTools

INPUT=$PROJ/text_data
OUTPUT=$PROJ/tappas_results

: <<'END'
#Standardized Euclidean Distance
echo "running SED"
python $SCRIPTS/standardized_euclidean_distance.py \
    --input $INPUT/transcript_tappas_amb_mean.tsv \
    --design $INPUT/df_genotypes.tsv \
    --ID Transcript \
    --figure $OUTPUT/SED_figure_mean_Amb.pdf \
    --SEDtoMean $OUTPUT/SED_samples_mean_Amb.tsv \
    --SEDpairwise $OUTPUT/SED_pairwise_mean_Amb.tsv \


END
#Log Transformation
echo "running log transformation"
python $SCRIPTS/log_and_glog_transformation.py \
    --input $INPUT/transcript_tappas_amb_mean.tsv \
    --design $INPUT/df_genotypes.tsv \
    --uniqID Transcript \
    --transformation 'log' \
    --log_base 'log' \
    --oname $OUTPUT/Log_mean_Amb.tsv \


usage: bland_altman_plot.py [-h] -i INPUT -d DESIGN -id UNIQID [-g GROUP] -f
                            BANAME -fd DISTNAME -fs FLAGSAMPLE -ff FLAGFEATURE
                            -pf PROPFEATURE -ps PROPSAMPLE
                            [-po PROCESSONLY [PROCESSONLY ...]]
                            [-rc RESIDCUTOFF] [-sfc SAMPLECUTOFF]
                            [-ffc FEATURECUTOFF] [--debug]
bland_altman_plot.py: error: the following arguments are required: -pf/--prop_feature, -ps/--prop_sample
#Bland Altman plots
echo "running BA plots on log transformed data"
python $SCRIPTS/bland_altman_plot.py \
    --input $OUTPUT/Log_mean_Amb.tsv \
    --design $INPUT/df_genotypes.tsv \
    --ID Transcript \
    --figure $OUTPUT/BA_figure_mean_Amb.pdf \
    --flag_dist $OUTPUT/BA_dist_flags_mean_Amb.pdf \
    --flag_sample $OUTPUT/BA_sample_flags_mean_Amb.tsv \
    --flag_feature $OUTPUT/BA_feature_flags_mean_Amb.tsv \
    --prop_feature $OUTPUT/BA_feature_prop_mean_Amb.tsv \
    --prop_sample $OUTPUT/BA_sample_prop_mean_Amb.tsv \
    --resid_cutoff 3 \
    --sample_flag_cutoff 0.2 \
    --feature_flag_cutoff 0.05

