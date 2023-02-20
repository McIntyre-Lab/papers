#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=prep
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/prep_rnaseq_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --array=47-56


date;hostname

### remove adapter sequences
### merge R1-R2 with bbmerge
### run fqSplitDups


## Using looping of the array task runs
## df (no header) has 336 rows
## want 6 rows per array task
        ## split into X array tasks:  336/X = 6 rows per array, X = 56

#Set the number of runs that each SLURM task should do
PER_TASK=6

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo -e "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    #Set directories
    PROJ=/ufrc/mcintyre/share/etoh_srna
    DATAIN=$PROJ/unzipped_rnaseq_FQ
    OUT=$PROJ/rnaseq

    LOGS=$OUT/dataPrepLogs
        mkdir -p $LOGS
    CUTADAPT=$OUT/cutadapt
        mkdir -p $CUTADAPT
    BBM=$OUT/cutadapt_bbmerge
        mkdir -p $BBM
    FQSPLIT=$OUT/cutadapt_bbmerge_fastqSplitDups
        mkdir -p $FQSPLIT

    ## DF:
    #        mel_R153_f_etoh_rep1,mel,R153,f,etoh,1,100,HLY23BBXX_UFE_172202,WF12,RAPiD-Genomics_HLY23BBXX_UFE_172202_P05_WF12_i5-509_i7-168_S456_L001

    ## Design file (for paired end reads)
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

    ## trim illumina truseq adapters
    ## running threaded, no quality trimming
    ## max error rate and min overlap -- parameters for identifying adapters, using trago/dros ones
	## specifying default values for max error rate and min overlap

    module purge
    module load cutadapt/2.1
    LOG=${LOGS}/${SAMPLEID}_${ID}_dataPrep.log
    ROZ=${LOGS}/roz_${SAMPLEID}_${ID}
        mkdir -p ${ROZ}

    echo -e "(1) STARTING CUTADAPT SAMPLE:  $SAMPLEID `date`\n" >${LOG}

    cutadapt \
        --cores=4 \
        -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;max_error_rate=0.1;min_overlap=3" \
        -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.1;min_overlap=3" \
        -o $CUTADAPT/${SAMPLEID}_${ID}_R1.fastq \
        -p $CUTADAPT/${SAMPLEID}_${ID}_R2.fastq \
        ${DATAIN}/${ID}_R1_*.fastq \
        ${DATAIN}/${ID}_R2_*.fastq \
        1>>${LOG}

    ## basic bbmerge to combine overlapping R1-R2 reads
    ## merges by overlap
	## pfilter=1 allows only perfect overlap
        ## minOverlap = 25% of RL

    echo -e "\nCALCULATE MIN OVERLAP FOR BBMERGE BASED ON READ LENGTH OF:    ${RL}  FOR SAMPLE:    ${SAMPLEID}" >> ${LOG}
    M=$(echo $RL*0.25 | bc)
    # MIN should be int
    MIN=${M%.*}
    echo -e "MIN OVERLAP IS:  ${MIN}\n" >>${LOG}
    echo -e "(2) STARTING BBMERGE FOR SAMPLE:    ${SAMPLEID} `date`\n" >>${LOG}

    module purge
    module load bbmap/38.44

    bbmerge.sh \
        in1=${CUTADAPT}/${SAMPLEID}_${ID}_R1.fastq \
        in2=${CUTADAPT}/${SAMPLEID}_${ID}_R2.fastq \
        out=${BBM}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq \
        outu1=${BBM}/${SAMPLEID}_${ID}_R1_unmerge_min${MIN}.fastq \
        outu2=${BBM}/${SAMPLEID}_${ID}_R2_unmerge_min${MIN}.fastq \
        outinsert=${LOGS}/${SAMPLEID}_${ID}_bbmerge_min${MIN}_read_sizes.txt \
        ihist=${LOGS}/${SAMPLEID}_${ID}_bbmerge_min${MIN}_ihist.txt \
        showhiststats=t \
        pfilter=1 \
        minOverlap=${MIN} \
        2>>${LOG}

    rm $CUTADAPT/${SAMPLEID}_${ID}_R1.fastq
    rm $CUTADAPT/${SAMPLEID}_${ID}_R2.fastq

    module purge
    module load python/2.7.14

    ## Drop reads that are less than 50% of starting readLength
    GE=${BBM}_GE_50perc_RL
        mkdir -p $GE
    LE=${BBM}_LESS_50perc_RL
        mkdir -p $LE
    WIDOW=${OUT}/widow_R1_R2_reads
        mkdir -p ${WIDOW}

    echo -e "\n(3) STARTING SPLIT OUT READS LESS 50% OF RL AND REMOVE EMPTY READ PAIRS FOR SAMPLE:    ${SAMPLEID} \n `date`" >>${LOG}
    python /ufrc/mcintyre/share/cegs2_MvsF_exp/scripts/bbM_50perc_readLen_and_drop_empty_reads_05amm.py \
        -i ${BBM}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq \
        -r1 ${BBM}/${SAMPLEID}_${ID}_R1_unmerge_min${MIN}.fastq \
        -r2 ${BBM}/${SAMPLEID}_${ID}_R2_unmerge_min${MIN}.fastq \
        -rl ${RL} \
        -g ${GE} \
        -l ${LE} \
        -o ${WIDOW}
        1>>${LOG}

    rm ${BBM}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq

    ## fqSplitdups on merged
    echo -e "\n(4) STARTING FQSQPLITDUPS ON BBMERGED - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${GE}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_bbmerge_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_${ID}_fqSplitDup_on_bbmerge_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_bbmerge_min${MIN}.log

    ## fqSplitDups SE on R1 widows
    echo -e "\n(5) STARTING FQSPLITDUPS ON R1 WIDOW - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    cp ${WIDOW}/${SAMPLEID}_${ID}_R1_unmerge_min${MIN}_widow.fastq \
        ${ROZ}/${SAMPLEID}_${ID}_R1_widow_min${MIN}.fastq

    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${ROZ}/${SAMPLEID}_${ID}_R1_widow_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_R1_widow_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_${ID}_fqSplitDup_on_R1_widow_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_R1_widow_min${MIN}.log

    ## fqSplitDups SE on R2 widows
    echo -e "\n(6) STARTING FQSPLITDUPS ON R2 WIDOW - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    cp ${WIDOW}/${SAMPLEID}_${ID}_R2_unmerge_min${MIN}_widow.fastq \
        ${ROZ}/${SAMPLEID}_${ID}_R2_widow_min${MIN}.fastq

    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${ROZ}/${SAMPLEID}_${ID}_R2_widow_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_R2_widow_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_${ID}_fqSplitDup_on_R2_widow_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_R2_widow_min${MIN}.log

    ## fqSplitDups PE on unmerged
    echo -e "\n(7) STARTING FQSPLITDUPS ON UNMERGED READS - PE SAMPLE:  ${SAMPLEID} `date`\n" >>${LOG}
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${WIDOW}/${SAMPLEID}_${ID}_R1_unmerge_min${MIN}_noEmpty.fastq \
        -r2 ${WIDOW}/${SAMPLEID}_${ID}_R2_unmerge_min${MIN}_noEmpty.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_unmerge_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_${ID}_fqSplitDup_on_unmerge_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_${ID}_fqSplitDup_on_unmerge_min${MIN}.log

    echo -e "\n(8) DATAPREP COMPLETED FOR SAMPLE: ${SAMPLEID} `date`" >>${LOG}

    rm ${BBM}/${SAMPLEID}_${ID}_R*_unmerge_min${MIN}.fastq
    rm ${GE}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq
    rm ${LE}/${SAMPLEID}_${ID}_bbmerge_min${MIN}.fastq

    rm -r ${ROZ}

done
