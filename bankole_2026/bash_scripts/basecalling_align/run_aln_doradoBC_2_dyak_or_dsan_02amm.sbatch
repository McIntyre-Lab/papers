#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=alnTO
#SBATCH --qos=mcintyre-cid
#SBATCH --account=mcintyre-cid
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/alnTO_%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=2-9

## align yak to yak
## san to san

## minimap2 fl + res  [ignoring unclass for now]


## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
PYCHOP_DOR=${PROJ}/dorado_basecalling/pychopper  ## decided not to mix guppy vs dorado basecalling
PYCHOP=$PROJ/pychopper

## create output dirs
ALN=${PROJ}/aln_dorReads_2_dyak2_or_dsan1
    mkdir -p $ALN

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
SAMPLE=${ARRAY[1]}
BC=${ARRAY[2]}
BCNUM=${ARRAY[3]}
TR=${ARRAY[4]}
TECH=${ARRAY[5]}
DATE=${ARRAY[6]}
ONT_DIR=${ARRAY[7]}

SAMPLEID=${SAMPLE}_TR${TR}

ROZ=$PROJ/roz_aln_${SAMPLEID}
    mkdir -p $ROZ


## cat together fl + rsc + unclass
cat ${PYCHOP}/fl_${SAMPLEID}.fq \
    ${PYCHOP}/rescue_${SAMPLEID}.fq \
    ${PYCHOP}/unclass_${SAMPLEID}.fq \
    > ${ROZ}/${SAMPLEID}.fastq


    module purge
    module load minimap/2.17 samtools/1.10  bedtools/2.29.2

    if [[ $SAMPLEID == dyak* ]] ; then
        REF_FA=/blue/mcintyre/share/references/dyak_Prin_Tai18E2_2.1/GCF_016746365.2/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna
        SAM=${ALN}/${SAMPLEID}_2_dyak2.sam
        BAM=${ROZ}/${SAMPLEID}_2_dyak2.bam
        BAMS=${ROZ}/${SAMPLEID}_2_dyak2_sort.bam
        TAG=dyak2
    else
        REF_FA=/blue/mcintyre/share/references/dsan_Prin_1.1/GCF_016746245.2/GCF_016746245.2_Prin_Dsan_1.1_genomic.fna
        SAM=${ALN}/${SAMPLEID}_2_dsan1.sam
        BAM=${ALN}/${SAMPLEID}_2_dsan1.bam
        BAMS=${ALN}/${SAMPLEID}_2_dsan1_sort.bam
        TAG=dsan1
    fi

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


