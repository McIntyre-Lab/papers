#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=prep
#SBATCH --account=secim
#SBATCH --qos=secim
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/prep2_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --array=49-120
#

date;hostname

### remove adapter sequences
### merge R1-R2 with bbmerge
### run fqSplitDups

## Using looping of the array task runs
## df (no header) has 120 rows
## want	5 rows per array task
        ## split into X array tasks:  120/X = 5 rows per array task, X=24


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

    #Set directories
    MCPYTHON=/ufrc/mcintyre/share/python.git
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    DATAIN=$PROJ/unzipped_FQ_RNAseq_novogene_1st_lane

    LOGS=$PROJ/dataPrepLogs_NG1
        mkdir -p $LOGS
    CUTADAPT=$PROJ/cutadapt_NG1
        mkdir -p $CUTADAPT
    BBM=$PROJ/cutadapt_bbmerge_NG1
        mkdir -p $BBM
    FQSPLIT=$PROJ/cutadapt_bbmerge_fastqSplitDups_NG1
        mkdir -p $FQSPLIT

    ## DF:  maize_rnaseq_design_file_noHeader.csv
	# sampleID,baseFQ,tubeID,Inbred_Line,Plant,Chamber,Ozone_Trt,PacBio
	# Mo17_P1_C1_Amb,zm_1,1,Mo17,P1,C1,Amb,
 	# Mo17_P2_C1_Amb,zm_2,2,Mo17,P2,C1,Amb,
	# Mo17_P3_C1_Amb,zm_3,3,Mo17,P3,C1,Amb,
	# Mo17_P4_C1_Amb,zm_4,4,Mo17,P4,C1,Amb,


    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_NG_1st_lane_w_sampleID_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    BASEFQ=${ARRAY[1]}
       ## = zm_113
 
    SAMPLEID=${SAMPLE}
    RL=150


    ## trim illumina truseq adapters
    ## running threaded, no quality trimming
    ## max error rate and min overlap -- parameters for identifying adapters, using trago/dros ones
	## specifying default values for max error rate and min overlap

    module purge
    module load cutadapt/2.1
    LOG=$LOGS/${SAMPLEID}_dataPrep.log

    echo -e "starting sample:  $SAMPLEID `date`\n" >>$LOG
    echo -e "starting cutadapt sample:  $SAMPLEID"

    cutadapt \
        --cores=8 \
        -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA;max_error_rate=0.1;min_overlap=3" \
        -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.1;min_overlap=3" \
        -o $CUTADAPT/${SAMPLEID}_R1.fastq \
        -p $CUTADAPT/${SAMPLEID}_R2.fastq \
        $DATAIN/${BASEFQ}_1.fq \
        $DATAIN/${BASEFQ}_2.fq \
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
        outinsert=$LOGS/${SAMPLEID}_bbmerge_read_sizes.txt \
        ihist=$LOGS/${SAMPLEID}_bbmerge_ihist.txt \
        showhiststats=t \
        pfilter=1 \
        minOverlap=${MIN} \
        2>>$LOG

    rm $CUTADAPT/${SAMPLEID}_R1.fastq
    rm $CUTADAPT/${SAMPLEID}_R2.fastq

    module purge
    module load python/2.7.14

    ## Drop reads that are less than 50% of starting readLength
    GE=$PROJ/cut_bbmerge_GE_50perc_RL 
        mkdir -p $GE
    LE=$PROJ/cut_bbmerge_LESS_50perc_RL
        mkdir -p $LE

    echo -e "starting split out reads less 50% of RL `date`" >>$LOG
    python $PROJ/scripts/bbM_50perc_read_length.py \
        -i $BBM/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        -r $RL \
        -g $GE \
        -l $LE \
        1>>$LOG
        
   rm $BBM/${SAMPLEID}_bbmerge_min${MIN}.fastq    
   
    ##fqSplitDups on merged
    echo -e "\n"
    echo -e "starting FqSplitDups on bbmerged - SE `date`\n" >> $LOG
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 $GE/${SAMPLEID}_bbmerge_min${MIN}.fastq \
        --outdir $FQSPLIT \
        -o $LOGS/${SAMPLEID}_fqSplitDup_on_bberge_min${MIN}_summary.csv \
        -t $FQSPLIT/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}_table.tsv \
        -g $LOGS/${SAMPLEID}_fqSplitDup_on_bbmerge_min${MIN}.log

    ## fqSplitDups on unmerged
    echo -e "\n"
    echo -e "starting FqSplitDups on unmerged reads - PE `date`\n" >>$LOG
    python $PROJ/scripts/fastqSplitDups_2MAI.py \
        -r1 $BBM/${SAMPLEID}_R1_unmerge_min${MIN}.fastq \
        -r2 $BBM/${SAMPLEID}_R2_unmerge_min${MIN}.fastq \
        --outdir $FQSPLIT \
        -o $LOGS/${SAMPLEID}_fqSplitDup_on_bbUnm_min${MIN}_summary.csv \
        -t $FQSPLIT/${SAMPLEID}_fqSplitDup_on_bbUnm_min${MIN}_table.tsv \
        -g $LOGS/${SAMPLEID}_fqSplitDupon_bbUnm_min${MIN}.log

    echo -e "completed for sample: $SAMPLEID `date`" >>$LOG

    rm $BBM/${SAMPLEID}_R*_unmerge_min${MIN}.fastq
    rm $GE/${SAMPLEID}_bbmerge_min${MIN}.fastq
    rm $LE/${SAMPLEID}_bbmerge_min${MIN}.fastq

done
