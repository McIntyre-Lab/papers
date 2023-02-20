#!/bin/sh


## Detection Above Background (DABG) for MACS2 peaks


module purge
module load python3/3.6


## Use previous temporary design file of non-input sample types (no replicates)
## Remakes it if it does not exist
DF=${ROZ}/df_Combine_Rep_Peaks.csv
if [ ! -e ${DF} ]; then
    awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
        cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}
fi

echo "
Detection above background... $(date)
"

## Begin summary counts file
SUMMARY=${DABG}/DABG_summary_counts.csv
echo "antibody,species,sex,treatment,consensus_peaks,called_all_reps,not_called_rep1(#DABG),not_called_rep2(#DABG),not_called_rep3(#DABG),detected_in_all" > ${SUMMARY}

## DABG for each sample	type
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...MACS2 ${TYPE}"

    ## Select all replicate count files
    ## Add a sampleID column to each and cat together
    REPS=$(ls ${CVR}/cvrg_cnts_${TYPE}_*_combined.csv)
    for file in ${REPS}; do
        NAME=$(basename ${file} _combined.csv | sed 's/cvrg_cnts_//')
        awk -v name=${NAME} '{if(NR==1){print "sampleID,"$0;} \
            else{print name","$0}}' $file > ${ROZ}/temp_${NAME}.csv
    done
    awk 'FNR==1 && NR!=1{next;}{print}' ${ROZ}/temp_* > ${ROZ}/combined_cnts.csv

    ## Get ASAW file with replicate information
    ASAW=${COMBINE}/nomod147ext_${TYPE}_asaw.csv

    ## Merge all coverage count files
    python3 ${SCRIPTS}/flag_DABG_04avn.py \
        -i ${ROZ}/combined_cnts.csv \
        -p 0.05 \
        -d ${ROZ}/${TYPE}.sql.db \
        -f ${ASAW} \
        -o ${DABG}/DABG_flag_merged_${TYPE}_coverage_cnts.tsv \
        -t ${DABG}/DABG_flag_merged_${TYPE}_coverage_cnts.full_table.tsv \
        >> ${SUMMARY} 2>${DABG}/DABG_flag_merged_${TYPE}_coverage_cnts.log

    ## Remove temporary coverage count files and sqlite3 database file
    rm ${ROZ}/temp_* ${ROZ}/combined_cnts.csv ${ROZ}/${TYPE}.sql.db
done

