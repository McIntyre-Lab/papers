#!/bin/bash

## Annotate MACS2 consensus peaks with genomic regions they overlap with

module purge
module load bedtools/2.28.0
module load python3/3.6

## Use previous temporary design file of non-input sample types (no replicates)
## Remakes it if it does not exist
DF=${ROZ}/df_Combine_Rep_Peaks.csv
if [ ! -e ${DF} ]; then
    awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
        cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}
fi

echo "
Annotate MACS2 consensus peaks... $(date)
"

## Annotate consensus peaks for each sample type
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...MACS2 ${TYPE}"

    ## Get consensus peaks BED file
    ASAW=${COMBINE}/nomod147ext_${TYPE}_asaw.bed

    ## Get fragments(exonic), intron, and intergenic regions
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        FRAGMENT=${FRAGMENT1}
        INTRON=${INTRON1}
        INTERGEN=${INTERGEN1}
    else
        FRAGMENT=${FRAGMENT2}
       	INTRON=${INTRON2}
       	INTERGEN=${INTERGEN2}
    fi


    ## Create genic bed file (all features other than intergenic combined)
    ## Check first if the final file exists already
    if [ ! -e ${FEATURES}/${SPECIES}_genic.merge.sorted.bed ]; then
        cat ${FEATURES}/${SPECIES}_*.unique.bed  ${FRAGMENT} ${INTRON} \
            | sort -k1,1 -k2,2n > ${ROZ}/${SPECIES}_genic.split.sorted.bed
        bedtools merge \
            -i ${ROZ}/${SPECIES}_genic.split.sorted.bed \
            > ${FEATURES}/${SPECIES}_genic.merge.sorted.bed
        ## Remove temporary file
        rm ${ROZ}/${SPECIES}_genic.split.sorted.bed
    fi

    ## Intersect consensus peaks with all feature types
    bedtools intersect -wo \
        -a ${ASAW} \
        -b ${FEATURES}/${SPECIES}_genic.merge.sorted.bed ${INTERGEN} \
        -names genic intergenic \
        > ${ANNOTATE}/${TYPE}_feature_intersect.tsv


    ## Get counts of genic, intergenic and border (both intergenic and genic intersections)
    VALUES=${ROZ}/${TYPE}_temp_values.tsv
    COLORS=${ROZ}/${TYPE}_temp_colors.tsv
    NAMES=${ROZ}/${TYPE}_temp_names.tsv
    python3 ${SCRIPTS}/count_peak_locations.py \
        -i ${ANNOTATE}/${TYPE}_feature_intersect.tsv \
        --output-counts ${VALUES} \
        --output-colors ${COLORS} \
        --output-names ${NAMES}

    ## Generate Pie Chart
    ## Possible values are genic, intergenic, or border
    echo "
Generating Pie Chart..."
    for FILETYPE in png tiff eps; do
        python3 ${SCRIPTS}/make_pie.py \
            -l $(cat ${NAMES}) \
            -v $(cat ${VALUES}) \
            -c $(cat ${COLORS}) \
            -o ${ANNOTATE}/${TYPE}_feature_intersect \
            -f ${FILETYPE}
    done

    ## Remove temporary files
    rm ${NAMES} ${COLORS} ${VALUES}

done
