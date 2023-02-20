#!/bin/sh


## Calculate fraction of reads in peaks (FRiP) for MACS2 peaks

module purge
module load python3/3.6
module load R/3.6

## Use previous temporary design file of non-input sample types (no replicates)
## Remakes it if it does not exist
DF=${ROZ}/df_Combine_Rep_Peaks.csv
if [ ! -e ${DF} ]; then
    awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
        cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}
fi

echo "
Fraction of reads in peaks... $(date)
"

## FRiP for each sample type
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...MACS2 ${TYPE}"

    ## Select all replicate count files and add a sampleID column
    REPS=$(ls ${CVR}/cvrg_cnts_${TYPE}_*_combined.csv)
    for file in ${REPS}; do
        NAME=$(basename ${file} _combined.csv | sed 's/cvrg_cnts_//')
        awk -v name=${NAME} '{if(NR==1){print "sampleID,"$0;} \
            else{print name","$0}}' ${file} > ${ROZ}/temp_${NAME}.csv
    done

    ## Select all corresponding input files and add sampleID column
    INPUTS=$(ls ${CVR}/${AB}_cvrg_cnts_I_${SPECIES}_${GENO}_${SEX}_${TRT}_*_combined.csv)
    for file in ${INPUTS}; do
        ## Keep the ${AB}_cvrg_cnts_ prefix so that the inputs can be distinguished later
        NAME=$(basename $file _combined.csv)
        awk -v name=${NAME} '{if(NR==1){print "sampleID,"$0;} \
            else{print name","$0}}' $file > ${ROZ}/temp_${NAME}.csv
    done

    ## Join together all coverage count files
    awk 'FNR==1 && NR!=1{next;}{print}' ${ROZ}/temp_*.csv > ${ROZ}/combined_cnts.csv

    ## Merge together all coverage count files and calculate FRiP
    python3 ${SCRIPTS}/calculate_frip.py \
        -i ${ROZ}/combined_cnts.csv \
        -o ${FRIP}/frip_${TYPE}.tsv \
        2>${FRIP}/frip_${TYPE}.log

    ## Remove temporary coverage count files
    rm ${ROZ}/temp_* ${ROZ}/combined_cnts.csv

done

## Summarize FRiP files
SUMMARY=${FRIP}/summary_frip_replicates.csv
awk 'FNR==1 && NR!=1 {next;}{print $1","$2","$3","$4","$5}' ${FRIP}/*.tsv \
    > ${SUMMARY}
