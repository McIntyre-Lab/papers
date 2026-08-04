#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=alnRMG
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/alnRMG_%A_%a.out
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
##SBATCH --array=8-13,20-25,32-37,44-49,56-61,68-73,80-85,92-97,104-109
#SBATCH --array=2

## aln mel to dmel8 and align sim samples to dsim202

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
PYCHOP=${PROJ}/dorado_basecalling/pychopper

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_bc_design_w_origDataPath_04amm.csv
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


## cat together fl + rsc + unclass
cat ${PYCHOP}/fl_${SAMPLEID}.fq \
    ${PYCHOP}/rescue_${SAMPLEID}.fq \
    ${PYCHOP}/unclass_${SAMPLEID}.fq \
    > ${ROZ}/${SAMPLEID}.fastq

    module purge
    module load minimap/2.17 samtools/1.10  bedtools/2.29.2

if [[ $SAMPLEID == mel* ]]; then
    REF_FA=/blue/mcintyre/share/references/dmel_fb650/dmel-all-chromosome-r6.50.fasta
    ALN=$PROJ/aln_mel_dorReads_2_dmel6
    SAM=${ALN}/${SAMPLEID}_2_dmel6.sam
    BAM=${ALN}/${SAMPLEID}_2_dmel6.bam
    BAMS=${ALN}/${SAMPLEID}_2_dmel6_sort.bam
    TAG=dmel6
else
    REF_FA=/blue/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
    ALN=$PROJ/aln_sim_dorReads_2_dsim2
    SAM=${ALN}/${SAMPLEID}_2_dsim2.sam
    BAM=${ROZ}/${SAMPLEID}_2_dsim2.bam
    BAMS=${ROZ}/${SAMPLEID}_2_dsim2_sort.bam
    TAG=dsim2
fi

    mkdir -p $ALN
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

