#!/bin/sh

## Call peaks with MACS2


module purge
module load gcc/5.2.0
module load macs/2.1.1


## Generate a temporary design file of all non-input samples
DF=${ROZ}/df_MACS2_Peak_Calling.csv
awk -F "," '$7 != "I"' ${DESIGN_FILE} > ${DF}


## Loop through samples in design file
echo "
Call MACS2 peaks for each sample... $(date)
"
for sample in $(cat ${DF}); do
    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $7}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')
    REP=$(echo ${sample} | awk -F "," '{print $5}')
    SAMPLEID=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}
    echo "    ...${SAMPLEID}"

    ## Get BAM files
    BAM=$BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam
    CONTROL=$BAMFILES/I_${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}_combined_bwa_PE_SE_uniq.sorted.bam

    ## MACS2 callpeaks
    macs2 callpeak \
        -t ${BAM} \
        -c ${CONTROL} \
        -q 0.05 \
        --name=nomod147ext_${SAMPLEID} \
        --outdir ${MACS} \
        --format=BAM \
        --gsize dm \
        --nomodel \
        --shift 73 \
        --extsize 147
done

