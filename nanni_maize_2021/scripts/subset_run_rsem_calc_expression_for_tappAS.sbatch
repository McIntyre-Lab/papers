#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=expMatrix
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/expMatrix_%A-%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --array=49-72

## Using looping of the array task runs
## df (no header) has 360 rows
## want 1 rows per array task
        ## split into X array tasks:  48/X = 1 rows per array task, X=48


#Set the number of runs that each SLURM task should do
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
    module load R rsem/1.2.31 star/2.7.0b

    ### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    REF=$PROJ/pacbio_analysis/subset_FSM_ISM_NIC_NNC/subset_references

    INPUT=$PROJ/cutadapt_bbmerge_fastqSplitDups
    INPUTNG=$PROJ/cutadapt_bbmerge_fastqSplitDups_NG1

    OUTDIR=$PROJ/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc
        mkdir -p $OUTDIR

    ## fq files in (1) cutadapt_bbmerge_fastqSplitDups AND (2) cutadapt_bbmerge_fastqSplitDups_NG1
    ## Hp301_P1_C1_Amb_techRep1_bbmerge_min37_distinct.fq  
    ## Hp301_P1_C1_Amb_techRep1_R1_unmerge_min37_distinct.fq  
    ## Hp301_P1_C1_Amb_techRep1_R2_unmerge_min37_distinct.fq


    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_rnaseq_samples_fix_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
	## sample = Mo17_P4_C5_Ele

    RSEMREF=$PROJ/pacbio_analysis/subset_FSM_ISM_NIC_NNC/subset_references/subset_fsm_ism_nic_nnc_rsem_index_STAR
    
    ROZ=$OUTDIR/roz_rsem_${SAMPLE}
        mkdir -p $ROZ


    ## (1) SE alignment of merged reads to reference PB filtered transcriptome ==> {SAMPLE}_techRep${TR}_bbmerge_min${MIN}_distinct.fq
    ## maize TRs are different lanes, RL is the same
    

    OUTSE=$OUTDIR/${SAMPLE}_SE

    cat ${INPUT}/${SAMPLE}_*_bbmerge_min*_distinct.fq ${INPUTNG}/${SAMPLE}_bbmerge_min*_distinct.fq > ${ROZ}/${SAMPLE}_bbmerge.fq
    READS=${ROZ}/${SAMPLE}_bbmerge.fq

    mkdir -p ${ROZ}_temp

    rsem-calculate-expression -p 4 \
        --temporary-folder ${ROZ}_temp \
        --keep-intermediate-files \
        --star \
        --output-genome-bam \
        --calc-ci \
        --estimate-rspd \
        ${READS} ${RSEMREF} ${OUTSE}

    ## (2) PE alignment of unmerged reads to reference PB filtered transcriptome ==> {SAMPLE}_techRep${TR}_R1/2_unmerge_min${MIN}_distinct.fq

    OUTPE=${OUTDIR}/${SAMPLE}_PE

    ## Concatenate reads and split empty mate pairs
    cat ${INPUT}/${SAMPLE}_*_R1_unmerge_min*_distinct.fq \
        ${INPUTNG}/${SAMPLE}_R1_unmerge_min*_distinct.fq \
        > ${ROZ}/${SAMPLE}_R1.fq

    cat ${INPUT}/${SAMPLE}_*_R2_unmerge_min*_distinct.fq \
        ${INPUTNG}/${SAMPLE}_R2_unmerge_min*_distinct.fq \
        > ${ROZ}/${SAMPLE}_R2.fq

    sh $PROJ/scripts/fastq_split_empty_mate_pairs.sh ${ROZ}/${SAMPLE}_R1.fq ${ROZ}/${SAMPLE}_R2.fq

    R1=${ROZ}/${SAMPLE}_R1.no_empty_pairs.fq
    R2=${ROZ}/${SAMPLE}_R2.no_empty_pairs.fq

    mkdir -p ${ROZ}_temp2

    rsem-calculate-expression -p 4 \
        --paired-end \
        --temporary-folder ${ROZ}_temp2 \
        --keep-intermediate-files \
        --star \
        --output-genome-bam \
        --calc-ci \
        --estimate-rspd \
        ${R1} ${R2} ${RSEMREF} ${OUTPE}

#    rm -r ${ROZ}

done
