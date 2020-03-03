#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=matrix
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/matrix_%A-%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --array=54,55,59,63,66,68,77,79,84,110
###SBATCH --array=1-4

## Using looping of the array task runs
## df (no header) has 360 rows
## want 5 rows per array task
        ## split into X array tasks:  120/X = 5 rows per array task, X=24


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
    module load R
    module load bowtie/0.12.9
    module load rsem/1.2.28

    ### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    REF=$PROJ/pacbio_analysis/references
    INPUT=$PROJ/cutadapt_bbmerge_fastqSplitDups

    RSEMREF=$REF/postFilter_PB_rsem_index_BT1

    OUTDIR=$PROJ/rsem_output_postFilter
        mkdir -p $OUTDIR
 
    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_rnaseq_samples_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    
    ROZ=$PROJ/rozRSE_${SAMPLE}_${RUN}
        mkdir -p $ROZ


    ## (1) SE alignment of merged reads to reference PB transcriptome ==> {SAMPLE}_techRep${TR}_bbmerge_min${MIN}_distinct.fq
    ## maize TRs are different lanes, RL is the same

    ## default aligner is bowtie1, or request --bowtie2  --star
    ## note jeremy used following options:
	##        --fragment-length-mean 80 \
	##        --fragment-length-sd 50 \

    OUTSE=$OUTDIR/${SAMPLE}_SE
        mkdir -p $OUTSE

    cat $INPUT/${SAMPLE}_*_bbmerge_min*_distinct.fq > $ROZ/${SAMPLE}_bbmerge.fq
    READS=$ROZ/${SAMPLE}_bbmerge.fq

    mkdir -p ${ROZ}_temp
    rsem-calculate-expression -p 4 \
        --temporary-folder ${ROZ}_temp \
        --calc-ci \
        --estimate-rspd \
        ${READS} ${RSEMREF} ${OUTSE}

    ## (2) PE alignment of unmerged reads to reference PB transcriptome ==>  {SAMPLE}_techRep${TR}_R1/2_unmerge_min${MIN}_distinct.fq

    OUTPE=$OUTDIR/${SAMPLE}_PE
        mkdir -p $OUTPE

    cat $INPUT/${SAMPLE}_*_R1_unmerge_min*_distinct.fq > $ROZ/${SAMPLE}_R1.fq
    cat $INPUT/${SAMPLE}_*_R2_unmerge_min*_distinct.fq > $ROZ/${SAMPLE}_R2.fq

    R1=$ROZ/${SAMPLE}_R1.fq
    R2=$ROZ/${SAMPLE}_R2.fq
    mkdir -p ${ROZ}_temp2
    rsem-calculate-expression -p 4 \
        --paired-end \
        --temporary-folder ${ROZ}_temp2 \
        --calc-ci \
        --estimate-rspd \
        ${R1} ${R2} ${RSEMREF} ${OUTPE}

    rm -r $ROZ

done
