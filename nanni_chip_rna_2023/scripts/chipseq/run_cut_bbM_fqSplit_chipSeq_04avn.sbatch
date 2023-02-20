#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=prep
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/prep_chip_%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-8


date;hostname

### remove adapter sequences
### merge R1-R2 with bbmerge
### run fqSplitDups


## Using looping of the array task runs
## df (no header) has 72 rows
## want 9 rows per array task
        ## split into X array tasks:  72/X = 9 rows per array, X = 8

#Set the number of runs that each SLURM task should do
PER_TASK=9

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
    DATAIN=$PROJ/unzipped_chipSeq_FQ
    OUT=$PROJ/chipseq

    LOGS=$PROJ/chipseq/dataPrepLogs
        mkdir -p $LOGS
    CUTADAPT=$PROJ/chipseq/cutadapt
        mkdir -p $CUTADAPT
    BBM=$PROJ/chipseq/cutadapt_bbmerge
        mkdir -p $BBM
    FQSPLIT=$PROJ/chipseq/cutadapt_bbmerge_fastqSplitDups
        mkdir -p $FQSPLIT

    ## DF:
    #      mel,r301,f,etoh,2,1a,I
    #      mel,r301,f,etoh,2,1a,K4
    #      mel,r301,f,etoh,2,1a,K27


    ## Design file (for paired end reads)
    DESIGN_FILE=$PROJ/design_files/design_etoh_srna_chipSeq_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SPECIES=${ARRAY[0]}
    GENO=${ARRAY[1]}
    SEX=${ARRAY[2]}
    TRT=${ARRAY[3]}
    REP=${ARRAY[4]}
    ID=${ARRAY[5]}
    AB=${ARRAY[6]}

    SAMPLEID=${AB}_${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}
    RL=100

    ## note that unzipped FQ files are prepended with Natalie's ID number

    ## trim illumina truseq adapters
    ## running threaded, no quality trimming
    ## max error rate and min overlap -- parameters for identifying adapters, using trago/dros ones
	## specifying default values for max error rate and min overlap

    module purge
    module load cutadapt/2.1
    LOG=$LOGS/${SAMPLEID}_dataPrep.log
    ROZ=${LOGS}/roz_${SAMPLEID}
        mkdir -p ${ROZ}

    echo -e "starting sample:  $SAMPLEID `date`\n" >>$LOG
    echo -e "starting cutadapt"

    cutadapt \
        --cores=8 \
        -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;max_error_rate=0.1;min_overlap=3" \
        -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.1;min_overlap=3" \
        -o $CUTADAPT/${SAMPLEID}_R1.fastq \
        -p $CUTADAPT/${SAMPLEID}_R2.fastq \
        $DATAIN/${ID}-${AB}_*_R1_*.fastq \
        $DATAIN/${ID}-${AB}_*_R2_*.fastq \
        1>>$LOG


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
        in1=$CUTADAPT/${SAMPLEID}_R1.fastq \
        in2=$CUTADAPT/${SAMPLEID}_R2.fastq \
        out=$BBM/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        outu1=$BBM/${SAMPLEID}_R1_unmerge_min${MIN}.fastq \
        outu2=$BBM/${SAMPLEID}_R2_unmerge_min${MIN}.fastq \
        outinsert=$LOGS/${SAMPLEID}_bbmerge_min${MIN}_read_sizes.txt \
        ihist=$LOGS/${SAMPLEID}_bbmerge_min${MIN}_ihist.txt \
        showhiststats=t \
        pfilter=1 \
        minOverlap=${MIN} \
        2>>$LOG

    rm $CUTADAPT/${SAMPLEID}_R1.fastq
    rm $CUTADAPT/${SAMPLEID}_R2.fastq

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
        -i ${BBM}/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        -r1 ${BBM}/${SAMPLEID}_R1_unmerge_min${MIN}.fastq \
        -r2 ${BBM}/${SAMPLEID}_R2_unmerge_min${MIN}.fastq \
        -rl ${RL} \
        -g ${GE} \
        -l ${LE} \
        -o ${WIDOW}
        1>>${LOG}

    rm $BBM/${SAMPLEID}_bbmerge_min${MIN}.fastq

    ## fqSplitdups on merged
    echo -e "\n(4) STARTING FQSQPLITDUPS ON BBMERGED - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${GE}/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}.log

    ## fqSplitDups on unmerged
    echo -e "\n(5) STARTING FQSPLITDUPS ON R1 WIDOW - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    cp ${WIDOW}/${SAMPLEID}_R1_unmerge_min${MIN}_widow.fastq \
        ${ROZ}/${SAMPLEID}_R1_widow_min${MIN}.fastq

    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${ROZ}/${SAMPLEID}_R1_widow_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_fqSplitDup_on_R1_widow_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_fqSplitDup_on_R1_widow_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_fqSplitDup_on_R1_widow_min${MIN}.log

    ## fqSplitDups SE on R2 widows
    echo -e "\n(6) STARTING FQSPLITDUPS ON R2 WIDOW - SE SAMPLE:  ${SAMPLEID} `date`\n" >> ${LOG}
    cp ${WIDOW}/${SAMPLEID}_R2_unmerge_min${MIN}_widow.fastq \
        ${ROZ}/${SAMPLEID}_R2_widow_min${MIN}.fastq

    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${ROZ}/${SAMPLEID}_R2_widow_min${MIN}.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_fqSplitDup_on_R2_widow_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_fqSplitDup_on_R2_widow_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_fqSplitDup_on_R2_widow_min${MIN}.log

    ## fqSplitDups PE on unmerged
    echo -e "\n(7) STARTING FQSPLITDUPS ON UNMERGED READS - PE SAMPLE:  ${SAMPLEID} `date`\n" >>${LOG}
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 ${WIDOW}/${SAMPLEID}_R1_unmerge_min${MIN}_noEmpty.fastq \
        -r2 ${WIDOW}/${SAMPLEID}_R2_unmerge_min${MIN}_noEmpty.fastq \
        --outdir ${FQSPLIT} \
        -o ${LOGS}/${SAMPLEID}_fqSplitDup_on_unmerge_min${MIN}_summary.csv \
        -t ${FQSPLIT}/${SAMPLEID}_fqSplitDup_on_unmerge_min${MIN}_table.tsv \
        -g ${LOGS}/${SAMPLEID}_fqSplitDup_on_unmerge_min${MIN}.log

    echo -e "\n(8) DATAPREP COMPLETED FOR SAMPLE: ${SAMPLEID} `date`" >>${LOG}

    rm ${BBM}/${SAMPLEID}_R*_unmerge_min${MIN}.fastq
    rm ${GE}/${SAMPLEID}_bbmerge_min${MIN}.fastq
    rm ${LE}/${SAMPLEID}_bbmerge_min${MIN}.fastq

    rm -r ${ROZ}

done
