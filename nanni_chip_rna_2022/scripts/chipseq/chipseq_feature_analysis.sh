#!/bin/sh


############################
###   Feature Analysis   ###
############################

## Examine ChIPseq reads over features for detection


## Set Output directories
export PILES=${ALN}/mpileups
    mkdir -p ${PILES}
export CVR=$PROJ/chipseq/coverage_cnts/features
     mkdir -p ${CVR}
export DABG=$PROJ/chipseq/detection_above_background/features
    mkdir -p ${DABG}
export FRIP=$PROJ/chipseq/frip/features
    mkdir -p ${FRIP}

## Get BAM file directory
export BAMFILES=${ALN}/bam_files

## Coverage counts on unique features
#source ${SCRIPTS}/coverage_counts_features_03avn.sh


## Fraction of reads in peaks (FRiP) of features
#source ${SCRIPTS}/frip_features_02avn.sh

## Check FRiP ratios
#sh ${SCRIPTS}/check_feature_frip_ratio.sh

## Plot FRiP values in box and whisker plot, split by ChIP marks
## Do for each species separately
module purge
module load R/3.6
#for SPECIES in ${SPECIES1} ${SPECIES2}; do
#    SUMMARY=${FRIP}/${SPECIES}_summary_frip_feature_replicates.csv
#    BOXPLOT=${FRIP}/${SPECIES}_summary_frip_feature
#    Rscript ${SCRIPTS}/boxplot_frip_features_02avn.R \
#        ${SUMMARY} ${BOXPLOT} ${SPECIES}
#done
Rscript ${SCRIPTS}/boxplot_frip_features_combine_species.R \
    ${FRIP}/${SPECIES1}_summary_frip_feature_replicates.csv \
    ${FRIP}/${SPECIES2}_summary_frip_feature_replicates.csv \
    ${FRIP}/${SPECIES1}_${SPECIES2}_summary_frip_feature \
    ${SPECIES1} ${SPECIES2}


## Detection above input of features
#source ${SCRIPTS}/detection_above_input_features_02avn.sh


## Average apnFF for detected feature (BED file for figures)
#source ${SCRIPTS}/avg_feature_apnFF_bed.sh


## Count detected features
#echo "Count Detected Features"
#source ${SCRIPTS}/count_detected_features_08avn.sh

## Make venn diagrams of counts
#source ${SCRIPTS}/featureType_venn_diagram_02avn.sh
