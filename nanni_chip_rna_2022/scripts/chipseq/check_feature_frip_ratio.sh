#!/bin/sh

## Check the relative proportions of K4 FRiP values between features

module purge
module load python/3.6

PROJ=/ufrc/mcintyre/share/etoh_srna
SCRIPTS=$PROJ/scripts/chipseq
FRIP=$PROJ/chipseq/frip/features
CHECK=$FRIP/check_frip_ratios
    mkdir -p ${CHECK}

## Check the proportions of K4 FRiP between features across each sample
for SPECIES in mel sim; do
    SUMMARY=${CHECK}/${SPECIES}_frip_ratios.csv
    STAT=${CHECK}/${SPECIES}_frip_ratios.stats.csv
    python ${SCRIPTS}/check_feature_frip_ratio.py \
        -i ${FRIP}/${SPECIES}_summary_frip_feature_replicates.csv \
        -o ${SUMMARY} \
        -s ${STAT}
done
