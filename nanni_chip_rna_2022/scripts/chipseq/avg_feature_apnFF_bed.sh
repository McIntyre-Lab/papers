#!/bin/sh


## Calculate average feature apnFF for detected features and output bed file
## Used in Jeremy's figure creation


module purge
module load python/3.6

echo "
Calculate average feature apnFF... $(date)
"

## Use previous temporary design file of non-input sample types (no replicates)
## Remakes it if it does not exist
DF=${ROZ}/df_Combine_Rep_Peaks.csv
if [ ! -e ${DF} ]; then
    awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
        cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}
fi

## Average feature apnFF for each sample type
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...features ${TYPE}"

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

    python ${SCRIPTS}/chip_avg_feature_apnFF_bed_02avn.py \
        -i ${DABG}/DAI_flag_merged_${TYPE}.csv \
        --fragment ${FRAGMENT} \
        --fusion ${FUSION} \
        --intron ${INTRON} \
        -d ${ROZ}/${TYPE}.avg.sql.db \
        -o ${DABG}/DAI_chip_${TYPE}_avg_apnFF.bed
    rm ${ROZ}/${TYPE}.avg.sql.db
done
