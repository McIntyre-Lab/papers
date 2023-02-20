#!/bin/sh


## Calculate fraction of reads in peaks (FRiP) for features

module purge
module load python/3.6


## Get feature lengths output file and add header
LENGTH=${FRIP}/total_lengths.csv
echo "species,feature_type,total_length" > ${LENGTH}

for SPECIES in ${SPECIES1} ${SPECIES2}; do

    ## Get correct fragment and intron files for species
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        FRAGMENT=${FRAGMENT1}
        INTRON=${INTRON1}
        FUSION=${FUSION1}
    else
        FRAGMENT=${FRAGMENT2}
        INTRON=${INTRON2}
        FUSION=${FUSION2}
    fi


    ## Pull out all count files for species
    ## Add a sampleID column to each and concatenate them
    for feature in 3UTR 5UTR fragments fusions intergenic introns TSS1kbWindow TSS300bpWindow; do
        SAMPLES=$(ls ${CVR}/${feature}/cvrg_cnts_*_${SPECIES}_*_combined.csv)
        for file in ${SAMPLES}; do
            NAME=$(basename ${file} _combined.csv | sed 's/cvrg_cnts_//')
            awk -v name=${NAME} '{if(NR==1){print "sampleID,"$0;} \
                else{print name","$0}}' $file > ${ROZ}/temp_${NAME}.csv
        done
        awk 'FNR==1 && NR!=1{next;}{print}' ${ROZ}/temp_*_${SPECIES}_*.csv > ${ROZ}/${SPECIES}_${feature}_combined_cnts.csv

        ## Calculate FRiP for each feature type
        echo "
FRiP for ${SPECIES} unique ${feature} coverage counts...
"
        python ${SCRIPTS}/calculate_frip.py \
            -i ${ROZ}/${SPECIES}_${feature}_combined_cnts.csv \
            -o ${ROZ}/temp_frip_${SPECIES}_${feature}.tsv \
            2>${FRIP}/frip_${SPECIES}_${feature}.log

        ## Add feature type to TSV
        awk -v feature=${feature} '{if(NR==1){print $0"\tfeature_type";} \
            else{print $0"\t"feature}}' ${ROZ}/temp_frip_${SPECIES}_${feature}.tsv \
            > ${FRIP}/frip_${SPECIES}_${feature}.tsv
    done

    ## Combine all FRiP values
    SUMMARY=${FRIP}/${SPECIES}_summary_frip_feature_replicates.csv
    awk 'FNR==1 && NR!=1 {next;}{print $1","$2","$3","$4","$5","$6}' \
        ${FRIP}/frip_${SPECIES}_*.tsv > ${SUMMARY}

    ## Normalize FRiP values by feature length
    ## Removed in this version of script -- not used

    ## Remove temporary coverage count and frip files
    rm ${ROZ}/temp_*.csv
    rm ${ROZ}/${SPECIES}_*_combined_cnts.csv
    rm ${ROZ}/temp_frip_${SPECIES}_*.tsv

done
