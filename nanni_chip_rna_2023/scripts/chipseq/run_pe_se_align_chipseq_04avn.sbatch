#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=bwa_PE_SE
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --mail-type=FAIL
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/bwa_PE_SE_chip.%A.%a.out
#SBATCH -t 3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=71

## Using looping of the array task runs
## df (no header) has 72 rows
## want 6 rows per array task
        ## split into X array tasks:  72/X = 6 rows per array, X = 12

#Set the number of runs that each SLURM task should do
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo -e "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    #Load modules
    module purge
    module load bwa/0.7.15 python3/3.6 samtools/1.9

    #Set directories
    PROJ=/ufrc/mcintyre/share/etoh_srna
    OUTPUT=$PROJ/chipseq/aln_genome
        mkdir -p $OUTPUT
    BAMFILES=$OUTPUT/bam_files
        mkdir -p $BAMFILES
    PARSE=$OUTPUT/bwa_parsed_sam_files
        mkdir -p $PARSE
    INPUT=$PROJ/chipseq/cutadapt_bbmerge_fastqSplitDups

    ROZ=$OUTPUT/temp_${SLURM_ARRAY_TASK_ID}
        mkdir -p $ROZ

    ## Design file
    DESIGN_FILE=$PROJ/design_files/design_etoh_srna_chipSeq_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SPECIES=${ARRAY[0]}
    GENO=${ARRAY[1]}
    SEX=${ARRAY[2]}
    TRT=${ARRAY[3]}
    REP=${ARRAY[4]}
    TUBE=${ARRAY[5]}
    AB=${ARRAY[6]}

    SAMPLEID=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}

    if [[ ${SPECIES} == "mel" ]]; then
        REF=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17_BWA
        FA=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
    else 
        REF=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02_BWA
        FA=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
    fi

    ## Cat together SE fq of the sample (merged)
    cat ${INPUT}/${SAMPLEID}_bbmerge_min*_*.fq \
        ${INPUT}/${SAMPLEID}_R?_widow_min*_*.fq > $ROZ/${SAMPLEID}_merged.fq

    ## Cat together PE fq of the sample (unmerged)
    cat ${INPUT}/${SAMPLEID}_R1_unmerge_min*_noEmpty_*.fq > $ROZ/${SAMPLEID}_R1.fq
    cat ${INPUT}/${SAMPLEID}_R2_unmerge_min*_noEmpty_*.fq > $ROZ/${SAMPLEID}_R2.fq

    ### (1) align
    ## SE alignment of bbmerge reads to reference genome
    bwa mem \
        -t 4 \
        -M ${REF} \
        ${ROZ}/${SAMPLEID}_merged.fq \
        > ${OUTPUT}/${SAMPLEID}_merged_bwa_SE.sam

    ### (2) align
    ### PE alignment of unmerged reads to reference genome
    bwa mem \
        -t 4 \
        -M ${REF} \
        $ROZ/${SAMPLEID}_R1.fq \
        $ROZ/${SAMPLEID}_R2.fq \
        > ${OUTPUT}/${SAMPLEID}_unmerged_bwa_PE.sam


    ### (3) Parse alignments
    python /ufrc/mcintyre/share/python.git/BWASplitSAM_08amp.py \
        -fq1 ${ROZ}/${SAMPLEID}_merged.fq \
        -s ${OUTPUT}/${SAMPLEID}_merged_bwa_SE.sam \
        --outdir ${PARSE}
    python /ufrc/mcintyre/share/python.git/BWASplitSAM_08amp.py \
        -fq1 $ROZ/${SAMPLEID}_R1.fq \
        -fq2 $ROZ/${SAMPLEID}_R2.fq \
        -s ${OUTPUT}/${SAMPLEID}_unmerged_bwa_PE.sam \
        --outdir ${PARSE}

    ### (4) Select uniq alignments
    ## sam to bam, sort and index the merge and unmerged 
    for i in merged unmerged; do
        if [[ ${i} == "merged" ]]; then
            aln="SE"
        else
            aln="PE"
        fi

        ## cat mapped and opposite sam files to get uniq mapping reads
        cat ${PARSE}/${SAMPLEID}_${i}_bwa_${aln}_mapped.sam \
            ${PARSE}/${SAMPLEID}_${i}_bwa_${aln}_opposite.sam \
            > ${PARSE}/${SAMPLEID}_${i}_bwa_${aln}_uniq.sam

        samtools view -ut ${FA}.fai -o $BAMFILES/${SAMPLEID}_${i}_bwa_${aln}_uniq.bam \
            ${PARSE}/${SAMPLEID}_${i}_bwa_${aln}_uniq.sam

        ## Remove full and parsed sam files
        rm ${OUTPUT}/${SAMPLEID}_${i}_bwa_${aln}.sam
        rm ${PARSE}/${SAMPLEID}_${i}_bwa_${aln}_*.sam
    done

    ### (5) Join PE and SE bam files and sort
    ## create single bam file for merged and unmerged reads, sort and index
    samtools merge -f \
        $BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.bam \
        $BAMFILES/${SAMPLEID}_unmerged_bwa_PE_uniq.bam \
        $BAMFILES/${SAMPLEID}_merged_bwa_SE_uniq.bam
    samtools sort -T $BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.tmp.sorted \
        -o $BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam \
        -O bam $BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.bam
    samtools index $BAMFILES/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam

    ## Remove unsorted bam files
    rm $BAMFILES/${SAMPLEID}_*_uniq.bam

    rm -r $ROZ
done
