#!/bin/sh

## Coverage Counts for chipseq

module purge
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
        FRAGMENT=${FRAGMENT1}
        INTRON=${INTRON1}
        FUSION=${FUSION1}
    else
        FASTA=${FASTA2}
        FRAGMENT=${FRAGMENT2}
        INTRON=${INTRON2}
        FUSION=${FUSION2}
    fi

    ## Generate mpileups if they do not already exist
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

    ## Get Coverage counts for alignments for each bed file
    for feature in fragments introns fusions 5UTR 3UTR TSS1kbWindow TSS300bpWindow intergenic; do
        if [[ ${feature} == "fragments" ]]; then
            BED=${FRAGMENT}
        elif [[ ${feature} == "introns" ]]; then
            BED=${INTRON}
        elif [[ ${feature} == "fusions" ]]; then
            BED=${FUSION}
        elif [[ ${feature} == "5UTR" ]]; then
            BED=${FEATURES}/${SPECIES}_5UTR.unique.bed
        elif [[ ${feature} == "3UTR" ]]; then
            BED=${FEATURES}/${SPECIES}_3UTR.unique.bed
        elif [[ ${feature} == "TSS1kbWindow" ]]; then
            BED=${FEATURES}/${SPECIES}_TSS1kbWindow.unique.bed
        elif [[ ${feature} == "TSS300bpWindow" ]]; then
            BED=${FEATURES}/${SPECIES}_TSS300bpWindow.unique.bed
        elif [[ ${feature} == "intergenic" ]]; then
            BED=${FEATURES}/${SPECIES}_intergenic.bed
        fi

        mkdir -p ${CVR}/${feature}/logs

        if [ ! -e ${CVR}/${feature}/cvrg_cnts_${SAMPLEID}_combined.csv ]; then
            python /ufrc/mcintyre/share/python.git/rpkm_calculate.py \
                -m ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup \
                -s ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam \
                -b ${BED} \
                -g ${CVR}/${feature}/logs/cvrg_cnts_${SAMPLEID}_combined.logfile \
                -o ${CVR}/${feature}/cvrg_cnts_${SAMPLEID}_combined.csv \
                -c
        fi
    done
    ## remove sorted sam file
    rm ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam

done

