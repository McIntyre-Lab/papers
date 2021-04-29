#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rsem_calc
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/rsem_calc_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --array=64,65,77

# Using looping of the array task runs
## df (no header) has 120 rows

#Set the number of runs that each SLURM task should do
##PER_TASK=5
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    echo "This is SLURM task $SLURM_ARRAY_TASK_ID, run number $RUN"

    module purge
    module load rsem/1.2.31 star/2.7.5c bwa/0.7.7 gcc/5.2.0 bowtie/0.12.9

    ### Set Directories
    PROJ=/blue/mcintyre/share/maize_ainsworth
    INPUT=$PROJ/cutadapt_bbmerge_fastqSplitDups
    INPUTNG=$PROJ/cutadapt_bbmerge_fastqSplitDups_NG1
    OUTDIR=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
    RSEMREF=${OUTDIR}/references/fsm_ism_nic_nnc_consol_junc_rsem_index_STAR

    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_rnaseq_samples_fix_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    
    ROZ=${OUTDIR}/roz_rsem_${SAMPLE}_${RUN}
        mkdir -p $ROZ

    ## (1) SE alignment of merged reads to transcriptome ==> {SAMPLE}_techRep${TR}_bbmerge_min${MIN}_distinct.fq

    ## default aligner is bowtie1, or request --bowtie2  --star

    ## Make single end read output directory
    OUTSE=${OUTDIR}/${SAMPLE}_SE
        mkdir -p ${OUTSE}

    ## Concatenate merged reads
    READS=${ROZ}/${SAMPLE}_merge.fq
    cat ${INPUT}/${SAMPLE}_*_bbmerge_min*_distinct.fq \
        ${INPUTNG}/${SAMPLE}_bbmerge_min*_distinct.fq \
        > ${READS}

    ## Set rsem temporary directory (is deleted at the end of rsem)
    mkdir -p ${ROZ}_temp

    ## Calculate expression for merged reads using STAR
    rsem-calculate-expression -p 4 \
        --temporary-folder ${ROZ}_temp \
        --keep-intermediate-files \
        --star \
        --output-genome-bam \
        --calc-ci \
        --estimate-rspd \
        ${READS} ${RSEMREF} ${OUTSE}

    ## (2) PE alignment of unmerged reads to reference PB transcriptome ==>  {SAMPLE}_techRep${TR}_R1/2_unmerge_min${MIN}_distinct.fq

    ## Make paired end read output directory
    OUTPE=${OUTDIR}/${SAMPLE}_PE
        mkdir -p ${OUTPE}

    ## Concatenate unmerged reads
    cat ${INPUT}/${SAMPLE}_*_R1_unmerge_min*_distinct.fq \
        ${INPUTNG}/${SAMPLE}_R1_unmerge_min*_distinct.fq \
        > ${ROZ}/${SAMPLE}_R1.fq
    cat ${INPUT}/${SAMPLE}_*_R2_unmerge_min*_distinct.fq \
        ${INPUTNG}/${SAMPLE}_R2_unmerge_min*_distinct.fq \
        > ${ROZ}/${SAMPLE}_R2.fq

    ## Split fastq with empty mate pairs
    sh $PROJ/scripts/fastq_split_empty_mate_pairs.sh ${ROZ}/${SAMPLE}_R1.fq ${ROZ}/${SAMPLE}_R2.fq
    R1=${ROZ}/${SAMPLE}_R1.no_empty_pairs.fq
    R2=${ROZ}/${SAMPLE}_R2.no_empty_pairs.fq


    ## Set rsem temporary directory (is deleted at the end of rsem)
    mkdir -p ${ROZ}_temp2

    ## Calculate expression for single end reads using STAR
    rsem-calculate-expression -p 4 \
        --paired-end \
        --temporary-folder ${ROZ}_temp2 \
        --keep-intermediate-files \
        --star \
        --output-genome-bam \
        --calc-ci \
        --estimate-rspd \
        ${R1} ${R2} ${RSEMREF} ${OUTPE}

    rm -r ${ROZ}

done
