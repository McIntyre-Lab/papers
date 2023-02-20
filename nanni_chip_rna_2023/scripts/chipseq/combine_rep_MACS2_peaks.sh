#!/bin/sh

## Combine all replicates (3 reps each) using average summit and average width (ASAW)

module purge
module load perl/5.24.1

## Generate temporary design file of non-input sample types (no replicates)
DF=${ROZ}/df_Combine_Rep_Peaks.csv
awk -F "," '$7 != "I"' ${DESIGN_FILE} | \
    cut -d "," -f 1,2,3,4,7 | sort | uniq > ${DF}

## Loop through sample types to combine replicate MACS2 peaks
echo "
Combine replicate peaks... $(date)
"
for sample in $(cat ${DF}); do
    ## Set sample variables and sample type
    AB=$(echo ${sample} | awk -F "," '{print $5}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')

    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...${TYPE}"

    PEAKS=$(ls ${MACS}/nomod147ext_${TYPE}_*_peaks.narrowPeak)

    ## Process each individual set of bed files first
    for INDIV in ${PEAKS}; do

        ## Get MACS2 peak bed files
        SUMMITS=${MACS}/$(basename ${INDIV} _peaks.narrowPeak)_summits.bed
        SINGLE=${MACS}/$(basename ${INDIV} _peaks.narrowPeak).sbed

        ## Combine summit bed files and narrow Peak bed file
        perl ${SCRIPTS}/macs_bed2sbed.pl \
            --bed ${INDIV} \
            --summit ${SUMMITS} \
            > ${SINGLE}
    done

    ## List all sbed files
    FILES=$(ls $MACS/nomod147ext_${TYPE}_*.sbed | \
        awk '{if(NR==1){all=$0;}else{all=all","$0;}}END{print all}')

    ## Combine peaks for all replicates with cmobine_peaks_v7.pl script
    ##      (v6 is listed for ChIP paper, but v7 is a more recent version)
    perl ${SCRIPTS}/combine_peaks_v7.pl \
        --asaw --files $FILES --reps 2 \
        > ${COMBINE}/nomod147ext_${TYPE}_asaw.csv

    ## Convert output csv to a bed file
    awk 'BEGIN{FS=","} {if(!/chrom/) print $2"\t"$3"\t"$4"\t"$1}' \
        ${COMBINE}/nomod147ext_${TYPE}_asaw.csv \
        > ${COMBINE}/nomod147ext_${TYPE}_asaw.bed

    ## Remove the single combined bed intermediate file (.sbed)
    rm ${SINGLE}
done
