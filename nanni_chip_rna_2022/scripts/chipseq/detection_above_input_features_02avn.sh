#!/bin/sh


## Detection Above Input for features


module purge
module load python/3.6

## Get median number of mapped reads across all samples in each species
## Use previously generated FRiP summary files (contains mapped reads)
MED1=$(awk -F "," '$1!="sampleID"{print $1"\t"$2}' \
    ${FRIP}/${SPECIES1}_summary_frip_feature_replicates.csv \
    | sort -k2,2n | uniq | \
    awk 'BEGIN{ \
            num=0; sum=0;} \
        {values[num++] = $2; sum+= $2;} \
        END{ \
            if(num%2==1){ \
                median=values[int(num/2)];} \
            else{ \
                median=(values[num/2]+values[(num/2)-1])/2} \
            printf "%.1f\n", median}')
echo "${SPECIES1} median mapped reads ff = ${MED1}
"

MED2=$(awk -F "," '$1!="sampleID"{print $1"\t"$2}' \
    ${FRIP}/${SPECIES2}_summary_frip_feature_replicates.csv \
    | sort -k2,2n | uniq | \
    awk 'BEGIN{ \
            num=0; sum=0;} \
        {values[num++] = $2; sum+= $2;} \
        END{ \
            if(num%2==1){ \
                median=values[int(num/2)];} \
            else{ \
                median=(values[num/2]+values[(num/2)-1])/2} \
            printf "%.1f\n", median}')
echo "${SPECIES2} median mapped reads ff = ${MED2}
"

## Use previous temporary design file of non-input sample types (no replicates)
## Remakes it if it does not exist
DF=${ROZ}/df_Combine_Rep_Peaks.csv
if [ ! -e ${DF} ]; then
    awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
        cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}
fi

echo "
Detection above input... $(date)
"

## Begin summary counts file
SUMMARY=${DABG}/DAI_summary_counts.csv
echo "sampleType,numFeatures,DAI_all,not_DAI_rep1,not_DAI_rep2,not_DAI_rep3" > ${SUMMARY}


## Detection above input for each sample type
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...features ${TYPE}"

    ## Combine all feature type coverage count files
    for file in ${CVR}/fragments/cvrg_cnts_${TYPE}_*_combined.csv; do
        NAME=$(basename ${file} _combined.csv | sed 's/cvrg_cnts_//')
        ALLCOMB=${ROZ}/all_features_${NAME}.csv
        HEADER=$(awk 'NR==1' ${file})
        for feature in 3UTR 5UTR fragments fusions intergenic introns TSS1kbWindow TSS300bpWindow; do
            awk 'NR!=1' ${CVR}/${feature}/cvrg_cnts_${NAME}_combined.csv \
                >> ${ALLCOMB}
        done
        sort -k1,1 -k2,2n ${ALLCOMB} >> ${ROZ}/all_features_${NAME}.sorted.csv
        echo ${HEADER} | cat - ${ROZ}/all_features_${NAME}.sorted.csv \
                > ${ROZ}/roz_file && mv ${ROZ}/roz_file ${ROZ}/all_features_${NAME}.sorted.csv
        rm ${ALLCOMB}
    done

    ## Combine Input features
    for file in ${CVR}/fragments/cvrg_cnts_I_${SPECIES}_${GENO}_${SEX}_${TRT}_*_combined.csv; do
        NAME=$(basename ${file} _combined.csv | sed 's/cvrg_cnts_//')
        ALLCOMB=${ROZ}/all_features_${NAME}.csv
        if [[ ! -e ${ROZ}/all_features_${NAME}.sorted.csv ]]; then
            HEADER=$(awk 'NR==1' ${file})
            for feature in 3UTR 5UTR fragments fusions intergenic introns TSS1kbWindow TSS300bpWindow; do
                awk 'NR!=1' ${CVR}/${feature}/cvrg_cnts_${NAME}_combined.csv \
                    >> ${ALLCOMB}
            done
            sort -k1,1 -k2,2n ${ALLCOMB} > ${ROZ}/all_features_${NAME}.sorted.csv
            echo ${HEADER} | cat - ${ROZ}/all_features_${NAME}.sorted.csv \
                > ${ROZ}/roz_file && mv ${ROZ}/roz_file ${ROZ}/all_features_${NAME}.sorted.csv 
            rm ${ALLCOMB}
        fi
    done

    ## Select all replicate count files
    ## Add a sampleID column to each and cat together
    REPS=$(ls ${ROZ}/all_features_${TYPE}_*.sorted.csv)
    for file in ${REPS}; do
        NAME=$(basename ${file} .sorted.csv | sed 's/all_features_//')
        awk -v name=${NAME} '{if(NR==1){print "sampleID,"$0;} \
            else{print name","$0}}' $file > ${ROZ}/temp_${NAME}.csv
    done

    ## Get proper median value for species
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        MED=${MED1}
    else
        MED=${MED2}
    fi    

    ## Flag detetection above input (DAI) for each rep
    for file in ${REPS}; do
        ## Get name of sample
        NAME=$(basename ${file} .sorted.csv | sed 's/all_features_//')

        ## Get corresponding input file
        INPUT=$(dirname ${file})/$(basename ${file} \
            | awk -F "_" 'BEGIN{OFS="_"}{$3 = "I"; print $0}')
        python3 ${SCRIPTS}/flag_apnFF_above_input.py \
            -c ${file} \
            -i ${INPUT} \
            -m ${MED} \
            -d ${ROZ}/temp_${NAME}.sql.db \
            -o ${ROZ}/flag_DAI_${NAME}.csv
        rm ${ROZ}/temp_${NAME}.sql.db
    done

    ## Flag feature presence in each feature type
    ## Feature must be DAI in at least 2 of the 3 replicates
    ## Merge all replicates for each sample type
    python3 ${SCRIPTS}/flag_DAI_presence.py \
        -n ${TYPE} \
        -1 ${ROZ}/flag_DAI_${TYPE}_rep1.csv \
        -2 ${ROZ}/flag_DAI_${TYPE}_rep2.csv \
        -3 ${ROZ}/flag_DAI_${TYPE}_rep3.csv \
        -d ${ROZ}/temp_${TYPE}.sql.db \
        -o ${DABG}/DAI_flag_merged_${TYPE}.csv \
        >> ${SUMMARY}

    ## Remove temporary database file
    rm ${ROZ}/temp_${TYPE}.sql.db
done

## Merge flags for all sample types
## Make stacked file with sample type name column
for SPECIES in ${SPECIES1} ${SPECIES2}; do
    FLAGS=$(ls ${DABG}/DAI_flag_merged_*_${SPECIES}_*.csv)
    for file in ${FLAGS}; do
        NAME=$(basename ${file} .csv | sed 's/DAI_flag_merged_//')
        awk -v name=${NAME} '{if(NR==1){print "sampleType,"$0;} \
            else{print name","$0}}' $file > ${ROZ}/temp_flags_${NAME}.csv
    done
    awk 'FNR==1 && NR!=1{next;}{print}' ${ROZ}/temp_flags_*.csv > ${ROZ}/stacked_flags.csv

    python ${SCRIPTS}/combine_DAI_flags.py \
        -i ${ROZ}/stacked_flags.csv \
        -o ${DABG}/${SPECIES}_DAI_summary_flags.csv
    
    ## Remove temporary files
    rm ${ROZ}/temp_flags_*.csv ${ROZ}/stacked_flags.csv
done
