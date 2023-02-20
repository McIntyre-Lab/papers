#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=combine_aln
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --mail-type=FAIL
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/combine_aln.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4gb
#SBATCH --array=1-47

## Using looping of the array task runs
## df (no header) has 48 rows
## want 1 rows per array task
        ## split into X array tasks:  48/X = 1 rows per array, X = 48

#Set the number of runs that each SLURM task should do
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo -e "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ )); do

    ## Load modules
    module purge
    module load bwa/0.7.7
    module load gcc/5.2.0
    module load python/2.7.10
    module load samtools/1.4
    module load bowtie/0.12.9

    ## Set directories
    PROJ=/ufrc/mcintyre/share/etoh_srna
    COMBINE=$PROJ/rnaseq/aln_genome/combine_techRep_bam_files
        mkdir -p ${COMBINE}
    BAMFILES=$PROJ/rnaseq/aln_genome/bam_files
    LOGS=${COMBINE}/logs
        mkdir -p $LOGS
    PILES=$PROJ/rnaseq/aln_genome/mpileups
        mkdir -p ${PILES}

    ## Design file
    DESIGN_FILE=$PROJ/design_files/design_ethanol_srna_RNAseq_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SPECIES=${ARRAY[1]}
    GENO=${ARRAY[2]}
    SEX=${ARRAY[3]}
    TRT=${ARRAY[4]}
    REP=${ARRAY[5]}
    WELLID=${ARRAY[7]}

    SAMPLEID=${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}

    if [[ ${SPECIES} == "mel" ]]; then
        FA=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
    else 
        FA=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
    fi

    ## Set log file for sample
    MYLOG=${LOGS}/${SAMPLEID}.log
    date
    echo "Combining techreps of sampleID ${SAMPLEID}

" >>${MYLOG}

##  Merge the SE and PE bam files for all techreps
    ## create single bam file for merged and unmerged reads, sort and index
    samtools merge -f \
        ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.bam \
        ${BAMFILES}/${SAMPLEID}_*_unmerge_*_bwa_PE_uniq.sorted.bam \
        ${BAMFILES}/${SAMPLEID}_*_widow_*_bwa_SE_uniq.sorted.bam \
        ${BAMFILES}/${SAMPLEID}_*_bbmerge_*_bwa_SE_uniq.sorted.bam
    samtools sort -T ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.tmp.sorted \
        -o ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam \
        -O bam ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.bam
    samtools index ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam

    ## Remove unsorted bam files
    rm ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.bam

    ## Create mpileup from merged bam files
    if [ ! -e ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup ]; then
        samtools mpileup \
            -d 1000000 \
            -f ${FA} \
            ${COMBINE}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam \
            > ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup
    fi

done
