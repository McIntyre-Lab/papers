#!/bin/bash

module purge
module load python/2.7

##########################################################
### The following environment variable need to be set ####
##########################################################
PROJ=${PROJ}                                          ####
OUTPUT=${OUTPUT}                                      ####
CVR=${CVR}                                            ####
ROZ=${ROZ}                                            ####
FEATURE=${FEATURE}                                    ####
SPECIES1=${SPECIES1}                                  ####
SPECIES2=${SPECIES2}                                  ####
UQFF=${UQFF}                                          ####
##########################################################

mkdir -p ${ROZ}

## Normalize APN by uq
## Average samples across all females and all males uq normalized APN
## Get ratio of female/male

for SPECIES in ${SPECIES1} ${SPECIES2}; do
    echo "
${SPECIES}...
"
    ## Normalize APN values with uq value
    python $PROJ/scripts/seq_xcrpt/apn_ff_wide_counts_03avn.py \
        -f ${UQFF} \
        -n uq_ff \
        -c ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_0_5.csv \
        -d ${ROZ}/${SPECIES}.ff.sql.db \
        -w ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_uq_ff.csv
    rm ${ROZ}/${SPECIES}.ff.sql.db

    ## Remove intermediate count file
    rm ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_0_5.csv

    ## Select fragments all on in one sex and all off in the other
    ## APN>5 flags
    awk -F "," 'NR==1 || ($26==1 && $29==1) || ($27==1 && $28==1)\
        {print $1","$26","$27","$28","$29}' \
        ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_uq_ff.csv \
        > ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_uq_ff.allsID_MF_apn5.csv

    ## Get ratios
    python $PROJ/scripts/seq_xcrpt/f_m_ratio_fragment_02avn.py \
        -i ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_uq_ff.csv \
        -o ${OUTPUT}/${SPECIES}_${FEATURE}_RNA

done

rm -r ${ROZ}
