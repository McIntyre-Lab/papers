#!/bin/sh

## Coverage Counts for chipseq

module purge
module load bwa/0.7.7
module load gcc/5.2.0
module load python/2.7.10
module load samtools/1.9


## Set mpileup directory
PILES=${ALN}/mpileups
    mkdir -p ${PILES}

## Create mpileup from merged bam files
echo "
Generate mpilups and coverage counts... $(date)
"
## Loop through each sample
for sample in $(cat ${DESIGN_FILE}); do

    ## Set sample variables and sampleID
    AB=$(echo ${sample} | awk -F "," '{print $7}')
    SPECIES=$(echo ${sample} | awk -F "," '{print $1}')
    GENO=$(echo ${sample} | awk -F "," '{print $2}')
    SEX=$(echo ${sample} | awk -F "," '{print $3}')
    TRT=$(echo ${sample} | awk -F "," '{print $4}')
    REP=$(echo ${sample} | awk -F "," '{print $5}')
    SAMPLEID=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}
    TYPE=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}
    echo "    ...${SAMPLEID}"

    ## Get reference FASTA file
    if [[ ${SPECIES} == ${SPECIES1} ]]; then
        FASTA=${FASTA1}
    else
        FASTA=${FASTA2}
    fi

    ## Generate mpilups if they do not already exist
    if [ ! -e ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup ]; then
        samtools mpileup \
            -d 1000000 \
            -f ${FASTA} \
            ${BAMFILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam \
            > ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup
    fi

    ## Convert sorted bam to sam for coverage counts if it does not already exist
    if [ ! -e ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam ]; then
        samtools view -h -o ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam \
            ${BAMFILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam
    fi

    ## Do coverage counts for consensus MACS2 peaks
    ## Check if sample is input, if so then there are two BED files
    if [[ ${AB} == "I" ]]; then
        BED=${COMBINE}/nomod147ext_*_${SPECIES}_${GENO}_${SEX}_${TRT}_asaw.bed
    else
        BED=${COMBINE}/nomod147ext_${TYPE}_asaw.bed
    fi

    for file in ${BED}; do
        ## Get prefix for input samples
        if [[ ${AB} == "I" ]]; then
            PREFIX=$(basename ${file} | awk -F "_" '{print $2}')_
        else
            PREFIX=""
        fi
        ## Get Coverage counts for alignments
        python /ufrc/mcintyre/share/python.git/tpm_calculate.py \
            -m ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup \
            -s ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam \
            -b ${file} \
            -g ${CVR}/logs/${PREFIX}cvrg_cnts_${SAMPLEID}_combined.logfile \
            -o ${CVR}/${PREFIX}cvrg_cnts_${SAMPLEID}_combined.csv \
            -c
    done

    ## remove sorted sam file
    rm ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam

done

