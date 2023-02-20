#!/bin/bash

module purge
module load python/3.6

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
##########################################################

mkdir -p ${ROZ}/counts

## (1) Convert coverage counts to wide format
## (2) Flag on/off for rnaseq

for SPECIES in ${SPECIES1} ${SPECIES2}; do
    echo "
${SPECIES}...
"
    ## (1) Get wide format counts
    ## Use python script to merge together all coverage count files
    cp ${CVR}/cvrg_cnts_${SPECIES}_*.csv ${ROZ}/counts/
    python $PROJ/scripts/merge_coverage_counts_apn.py \
        -d ${ROZ}/counts \
        -o ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.csv

    ## Clear temporary counts directory
    rm ${ROZ}/counts/*

    ## (2) Flag on/off for rnaseq

    ## Make design file
    DESIGN=${ROZ}/${SPECIES}.df.csv
    echo "sampleID,sampleType,species_sex,species,genotype,sex,treatment,rep" > ${DESIGN}
    awk -F "," -v species=${SPECIES} '$2==species{OFS=","; print $2"_"$3"_"$4"_"$5"_rep"$6,\
        $2"_"$3"_"$4"_"$5,$2"_"$4,$2,$3,$4,$5,$6}' $PROJ/design_files/design_ethanol_srna_RNAseq.csv \
        >> ${DESIGN}

    ## Flag on/off apn>0 and apn>5
    for apn in 0 5; do
        python $PROJ/scripts/flag_event_detection_03avn.py \
            -i ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.csv \
            -d ${DESIGN} \
            -g sampleType \
            -a ${apn} -p 0.5 \
            -f species_sex \
            -s species_sex \
            -o ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.on_off_${apn}.csv

    ## (3) Count flags and select features all on in one sex and all off in the other
        ## Count on/off m/f crosstabs
        python $PROJ/scripts/count_on_off_flag.py \
            -f ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.on_off_${apn}.csv \
            -o ${OUTPUT}/merged_${SPECIES}_${FEATURE}_coverage_cnts.on_off_${apn}.counts.txt
    done

    ## Merge flags and counts
    python $PROJ/scripts/merge_on_off_flag.py \
        -1 ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.on_off_0.csv \
        -2 ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.on_off_5.csv \
        -c ${ROZ}/merged_${SPECIES}_${FEATURE}_coverage_cnts.csv \
        -d ${ROZ}/${SPECIES}.merge_0_5.sql.db \
        -o ${OUTPUT}/${SPECIES}_${FEATURE}_on_off_0_5.csv
    rm ${ROZ}/${SPECIES}.merge_0_5.sql.db

done

rm -r ${ROZ}
