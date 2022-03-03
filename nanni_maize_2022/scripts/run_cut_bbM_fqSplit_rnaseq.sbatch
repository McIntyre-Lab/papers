#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=prep
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/prep_%A.%a.out
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-18
#

date;hostname

### remove adapter sequences
### merge R1-R2 with bbmerge
### run fqSplitDups

## Using looping of the array task runs
## df (no header) has 360 rows
## want	20 rows per array task
        ## split into X array tasks:  360/X = 20 rows per array task, X=18


#Set the number of runs that each SLURM task should do
PER_TASK=20


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
    MCPYTHON=/ufrc/mcintyre/share/python.git
    PROJ=/ufrc/mcintyre/share/etoh_srna
    DATAIN=$PROJ/original_data/RNAseq/Plate3-4

    LOGS=$PROJ/rnaseq/dataPrepLogs
        mkdir -p $LOGS
    CUTADAPT=$PROJ/rnaseq/cutadapt
        mkdir -p $CUTADAPT
    BBM=$PROJ/rnaseq/cutadapt_bbmerge
        mkdir -p $BBM
    FQSPLIT=$PROJ/rnaseq/cutadapt_bbmerge_fastqSplitDups
        mkdir $FQSPLIT

    ## DF:  maize_rnaseq_design_file_noHeader.csv
    #	sampleID,baseFQ,RAPiD_Genomics_Sample_Code,tubeID,Inbred_Line,Plant,Chamber,Ozone_Trt,plate,techRep
    #   Mo17_P1_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC12_i5-501_i7-43_S567_L003,UFE_032201_P003_WC12,1,Mo17,P1,C1,Amb,1,1
    #   Mo17_P1_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC12_i5-501_i7-43_S687_L004,UFE_032201_P003_WC12,1,Mo17,P1,C1,Amb,1,2
    #   Mo17_P1_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC12_i5-501_i7-43_S807_L005,UFE_032201_P003_WC12,1,Mo17,P1,C1,Amb,1,3
    #   Mo17_P2_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC07_i5-501_i7-62_S562_L003,UFE_032201_P003_WC07,2,Mo17,P2,C1,Amb,1,1
    #   Mo17_P2_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC07_i5-501_i7-62_S682_L004,UFE_032201_P003_WC07,2,Mo17,P2,C1,Amb,1,2
    #   Mo17_P2_C1_Amb,RAPiD-Genomics_F086_UFE_032201_P003_WC07_i5-501_i7-62_S802_L005,UFE_032201_P003_WC07,2,Mo17,P2,C1,Amb,1,3


    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_renaseq_design_file_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    BASEFQ=${ARRAY[1]}
       ## = RAPiD-Genomics_F086_UFE_032201_P003_WA01_i5-501_i7-59_S532_L003
    TR=${ARRAY[9]}
 
    SAMPLEID=${SAMPLE}_techRep${TR}
    RL=150


    ## trim illumina truseq adapters
    ## running threaded, no quality trimming
    ## max error rate and min overlap -- parameters for identifying adapters, using trago/dros ones
	## specifying default values for max error rate and min overlap

    module purge
    module load cutadapt/2.1
    LOG=$LOGS/${SAMPLEID}_dataPrep.log

    echo -e "starting sample:  $SAMPLEID `date`\n" >>$LOG
    echo -e "starting cutadapt"

    cutadapt \
        --cores=8 \
        -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;max_error_rate=0.1;min_overlap=3" \
        -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.1;min_overlap=3" \
        -o $CUTADAPT/${SAMPLEID}_R1.fastq \
        -p $CUTADAPT/${SAMPLEID}_R2.fastq \
        $DATAIN/${BASEFQ}_R1_*.fastq.gz \
        $DATAIN/${BASEFQ}_R2_*.fastq.gz \
        1>>$LOG


    ## basic bbmerge to combine overlapping R1-R2 reads
    ## merges by overlap
	## pfilter=1 allows only perfect overlap
        ## minOverlap = 25% of RL

    echo -e "\n"
    echo -e "calculate minOverlap for bbmerge based on readLength of: $RL" >> $LOG
    M=$(echo $RL*0.25 | bc)
    # MIN should be int
    MIN=${M%.*}
    echo -e "minOverlap is:  $MIN\n" >>$LOG
    echo -e "starting bbmerge `date`\n" >>$LOG

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


    ### TO ADD:  drop reads that are less than 50% of starting readLength






    module purge
    module load python/2.7.6

    ## fqSplitdups on merged
    echo -e "\n"
    echo -e "starting FqSplitDups on bbmerged reads - SE `date\n`" >>$LOG
    $MCPYTHON/fastqSplitDups.py \
        -r1 $BBM/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        --outdir $FQSPLIT \
        -o $LOGS/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}_summary.csv \
        -t $FQSPLIT/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}_table.tsv \
        -g $LOGS/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}.log

    ## fqSplitDups on unmerged
    echo -e "\n"
    echo -e "starting FqSplitDups on unmerged reads - PE `date\n`" >>$LOG
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 $BBM/${SAMPLEID}_R1_unmerge_min${MIN}.fastq \
        -r2 $BBM/${SAMPLEID}_R2_unmerge_min${MIN}.fastq \
        --outdir $FQSPLIT \
        -o $LOGS/${SAMPLEID}_fqSplitDup_on_bbUnm_min${MIN}_summary.csv \
        -t $FQSPLIT/${SAMPLEID}_fqSplitDup_on_bbUnm_min${MIN}_table.tsv \
        -g $LOGS/${SAMPLEID}_fqSplitDupon_bbUnm_min${MIN}.log

    echo -e "completed for sample: $SAMPLEID `date`" >>$LOG

done
