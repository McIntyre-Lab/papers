#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=aln
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --mail-type=FAIL
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/aln_jxn_and_genome.%A.%a.out
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-84

## Using looping of the array task runs
## df (no header) has 336 rows
## want 4 rows per array task
        ## split into X array tasks:  336/X = 4 rows per array, X = 84

#Set the number of runs that each SLURM task should do
PER_TASK=4

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
    module load samtools/1.4
    module load bowtie/0.12.9

    ## Set directories
    PROJ=/ufrc/mcintyre/share/etoh_srna
    GENOMEOUT=${PROJ}/rnaseq/aln_genome
        mkdir -p ${GENOMEOUT}
    JXNOUT=${PROJ}/rnaseq/aln_junction
        mkdir -p ${JXNOUT}
    GENOMEBAM=${GENOMEOUT}/bam_files
        mkdir -p ${GENOMEBAM}
    JXNBAM=${JXNOUT}/bam_files
        mkdir -p ${JXNBAM}
    GENOMEPARSE=${GENOMEOUT}/bwa_parsed_sam_files
        mkdir -p ${GENOMEPARSE}
    GENOMELOGS=${GENOMEOUT}/logs
        mkdir -p ${GENOMELOGS}
    JXNLOGS=${JXNOUT}/logs
        mkdir -p ${JXNLOGS}
    INPUT=$PROJ/rnaseq/cutadapt_bbmerge_fastqSplitDups

    ## Design file
    DESIGN_FILE=$PROJ/design_files/design_PBrna_sampleID_readLength_fq_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SPECIES=${ARRAY[0]}
    GENO=${ARRAY[1]}
    SEX=${ARRAY[2]}
    TRT=${ARRAY[3]}
    REP=${ARRAY[4]}
    RL=${ARRAY[5]}
    RUNID=${ARRAY[6]}
    WELL=${ARRAY[7]}
    ID=${ARRAY[8]}

    SAMPLEID=${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}

    TMPDIR=${PROJ}/rnaseq/roz_aln_${SAMPLEID}_${ID}
        mkdir -p ${TMPDIR}

    ## Get value of min read threshold from read prep (25% RL)
    MIN=$(echo "${RL}/4" | bc)

    if [[ ${SPECIES} == "mel" ]]; then
        BWAREF=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17_BWA
        JXNFA=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_junctions_${RL}bp.fa
        FA=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
        REF=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_junctions_${RL}bp_BT1
    else 
        BWAREF=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02_BWA
        FA=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
        if [[ ${RL} == "150" ]]; then
            REF=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_junctions_BT1
            JXNFA=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_junctions.fa
        else
            REF=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_junctions_100bp_BT1
            JXNFA=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_junctions_100bp.fa
        fi
    fi

    ## Set junctions log file for sample
    ## must set ALNLOGS for alignment_functions.sh script
    ALNLOGS=${JXNLOGS}
    MYLOG=${JXNLOGS}/${SAMPLEID}_${ID}.log
    ALN_ERROR_LOG=${JXNLOGS}/size_errors.txt
    date
    echo "sampleID ${SAMPLEID}_${ID} with read length ${RL}

" >${MYLOG}

## (1) Align reads SE to junctions with bowtie1
    ## Store number of processors that will be used, i.e. the number of files we will split into to run LAST on
    NUMPROCS=1
    for i in bbmerge R1_unmerge R2_unmerge R1_widow R2_widow; do
        if [[ ${i} == "R1_unmerge" || ${i} == "R2_unmerge" ]]; then
            NAME=${SAMPLEID}_${ID}_${i}_min${MIN}_noEmpty
        else
            NAME=${SAMPLEID}_${ID}_${i}_min${MIN}
        fi

        ## Set reads
        READS=${INPUT}/${NAME}_distinct.fq

        echo "
Align ${i} reads SE to junctions with bowtie 1
        Junctions: ${REF}
        Reads: ${READS}
" >> ${MYLOG}
        source $PROJ/scripts/alignment_functions.sh
        bowtie_se_uniq

        ## Go up one directory (out of $TMPDIR)
        ## alignment_functions.sh has cd $TMPDIR
        cd ..

        if [ ! -e ${JXNFA}.fai ]; then
            samtools faidx ${JXNFA}
        fi
        samtools view -ut ${JXNFA}.fai -o ${TMPDIR}/${NAME}_aln_bt_junctions.bam \
            ${TMPDIR}/${NAME}_aln_bt.sam
        echo "samtools sort junctions"
        samtools sort -T ${TMPDIR}/${NAME}_aln_bt_junctions.tmp.sorted \
            -o ${JXNBAM}/${NAME}_aln_bt_junctions.sorted.bam \
            -O bam ${TMPDIR}/${NAME}_aln_bt_junctions.bam
        echo "samtools index junctions"
        samtools index ${JXNBAM}/${NAME}_aln_bt_junctions.sorted.bam

        echo "BT1 to junctions complete for ${NAME}" >>${MYLOG}

    done

## (2) Align all reads to genome
## merged and widow reads aligned SE, unmerged reads aligned PE

## Set genome log file for sample
    MYLOG=${GENOMELOGS}/${SAMPLEID}_${ID}.log
    date
    echo "sampleID ${SAMPLEID}_${ID} with read length ${RL}

" >${MYLOG}

    for i in bbmerge widow unmerge; do
        ## Concatenate reads
        NAME=${SAMPLEID}_${ID}_${i}_min${MIN}
        if [[ ${i} == "bbmerge" ]]; then
            aln="SE"
            READS=${INPUT}/${NAME}_distinct.fq
        elif [[ ${i} == "widow" ]]; then
            aln="SE"
            READS=${TMPDIR}/${NAME}_distinct.fq
            cat ${INPUT}/${SAMPLEID}_${ID}_R?_widow_min${MIN}_distinct.fq \
                > ${READS}
        else
            aln="PE"
            READS=$(ls ${INPUT}/${SAMPLEID}_${ID}_R?_unmerge_min${MIN}_noEmpty_distinct.fq)
        fi
        echo "
Align ${i} reads ${aln} to genome with bwa-mem
        Reference: ${BWAREF}
" >> ${MYLOG}
        bwa mem \
            -t 4 \
            -M ${BWAREF} \
            ${READS} \
            > ${GENOMEOUT}/${NAME}_bwa_${aln}.sam 2>>${MYLOG}
    done

## (3) Parse alignments
module load python3/3.6
    python3 /ufrc/mcintyre/share/python.git/BWASplitSAM_09avn.py \
        -fq1 ${INPUT}/${SAMPLEID}_${ID}_bbmerge_min${MIN}_distinct.fq \
        -s ${GENOMEOUT}/${SAMPLEID}_${ID}_bbmerge_min${MIN}_bwa_SE.sam \
        --outdir ${GENOMEPARSE} 2>>${MYLOG}
    python3 /ufrc/mcintyre/share/python.git/BWASplitSAM_09avn.py \
        -fq1 ${TMPDIR}/${SAMPLEID}_${ID}_widow_min${MIN}_distinct.fq \
        -s ${GENOMEOUT}/${SAMPLEID}_${ID}_widow_min${MIN}_bwa_SE.sam \
        --outdir ${GENOMEPARSE} 2>>${MYLOG}
    python3 /ufrc/mcintyre/share/python.git/BWASplitSAM_09avn.py \
        -fq1 ${INPUT}/${SAMPLEID}_${ID}_R1_unmerge_min${MIN}_noEmpty_distinct.fq \
        -fq2 ${INPUT}/${SAMPLEID}_${ID}_R2_unmerge_min${MIN}_noEmpty_distinct.fq \
        -s ${GENOMEOUT}/${SAMPLEID}_${ID}_unmerge_min${MIN}_bwa_PE.sam \
        --outdir ${GENOMEPARSE} 2>>${MYLOG}


## (4) Get uniquely mapping reads for genome alignments (mapped + opposite)
    ## sam to bam, sort and index the merged and unmerged 
    for i in bbmerge widow unmerge; do

        NAME=${SAMPLEID}_${ID}_${i}_min${MIN}

        if [[ ${i} == "unmerge" ]]; then
            aln="PE"
        else
            aln="SE"
        fi

        ## cat mapped and oppposite sam files to get uniq mapping reads
        cat ${GENOMEPARSE}/${NAME}_bwa_${aln}_mapped.sam \
            ${GENOMEPARSE}/${NAME}_bwa_${aln}_opposite.sam \
            > ${GENOMEPARSE}/${NAME}_bwa_${aln}_uniq.sam

        echo "samtools view genome"
        samtools view -ut ${FA}.fai -o ${GENOMEBAM}/${NAME}_bwa_${aln}_uniq.bam \
            ${GENOMEPARSE}/${NAME}_bwa_${aln}_uniq.sam
        echo "samtools sort genome"
        samtools sort -T ${GENOMEBAM}/${NAME}_bwa_${aln}_uniq.tmp.sorted \
            -o ${GENOMEBAM}/${NAME}_bwa_${aln}_uniq.sorted.bam \
            -O bam ${GENOMEBAM}/${NAME}_bwa_${aln}_uniq.bam
        echo "samtools index genome"
        samtools index ${GENOMEBAM}/${NAME}_bwa_${aln}_uniq.sorted.bam

       	## Remove the full unparsed sam file
        rm ${GENOMEOUT}/${NAME}_bwa_${aln}.sam
    done

    ## Remove unsorted bam files
    rm ${GENOMEBAM}/${SAMPLEID}_${ID}_*_uniq.bam

    rm -r ${TMPDIR}
done
