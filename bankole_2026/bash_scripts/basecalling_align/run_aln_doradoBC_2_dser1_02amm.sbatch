#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=alnK
#SBATCH --qos=mcintyre-cid
#SBATCH --account=mcintyre-cid
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/alnK_%A_%a.out
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=2-49
##SBATCH --array=2


## align to serrata genome


## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
PYCHOP=${PROJ}/dorado_basecalling/pychopper

## create output dirs
ALN=${PROJ}/aln_ser_dorReads_2_dser1
    mkdir -p $ALN

REF_FA=/blue/mcintyre/share/references/dser1.1/GCF_002093755.2/GCF_002093755.2_Dser1.1_genomic.fna

## Get design file information (2 - 49)
DESIGN_FILE=$PROJ/design_files/design_file_sample2BC_w_origDataPath_04amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

ROZ=$PROJ/roz_aln_${SAMPLEID}
    mkdir -p $ROZ

## cat together fl + rsc + unlcass
cat ${PYCHOP}/fl_${SAMPLEID}.fq \
    ${PYCHOP}/rescue_${SAMPLEID}.fq \
    ${PYCHOP}/unclass_${SAMPLEID}.fq \
    > ${ROZ}/${SAMPLEID}.fastq


    module purge
    module load minimap/2.17 samtools/1.10  bedtools/2.29.2

    SAM=${ALN}/${SAMPLEID}_2_dser1.sam

## map to ref
    minimap2 \
        -t 12 \
        -a \
        -x splice \
        --secondary=yes \
        -N 200 \
        -p 0.9 \
        -C 5 \
	${REF_FA} \
        ${ROZ}/${SAMPLEID}.fastq \
        > ${SAM}


module purge
module load samtools/1.10

TAG=dser1
## remove reads with mapping quality of 0
    samtools view \
        -h \
        -q 1 \
	${SAM} \
        > ${ALN}/${SAMPLEID}_2_${TAG}_noMQ0.sam

## exclude header, include reads without flag 2052 (omit unmapped and supplementary)
    samtools view \
        -h \
        -F 2052 \
        ${ALN}/${SAMPLEID}_2_${TAG}_noMQ0.sam \
        > ${ALN}/${SAMPLEID}_2_${TAG}_noMQ0_uniq.sam

## count mapped reads
    MAP=$(samtools view -c ${ALN}/${SAMPLEID}_2_${TAG}_noMQ0_uniq.sam)
    CNT=$(awk 'END {print NR/4}' ${ROZ}/${SAMPLEID}.fastq)

    echo "sample,input_read_num,mapped_read_num,prop_mapped" \
        > ${ALN}/${SAMPLEID}_2_${TAG}_mapped_read_cnts.csv
    echo "${SAMPLEID},${CNT},${MAP},$(echo "${MAP} ${CNT}" | awk '{print $1/$2}')" \
        >> ${ALN}/${SAMPLEID}_2_${TAG}_mapped_read_cnts.csv


## removed unneeded files
    rm $SAM
    rm ${ALN}/${SAMPLEID}_2_${TAG}_noMQ0.sam
    rm -r $ROZ


