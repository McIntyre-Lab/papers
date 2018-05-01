#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=aln_splicing
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/event_analysis/scripts/SLURM_LOGS/out.aln_splicing.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6gb
#SBATCH --array=1-4
#


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}


module load bowtie/0.12.9
module load python/2.7.6

## Store number of processors that will be used, i.e. the number of files we will split into to run LAST on
    NUMPROCS=4

## Set Directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check/event_analysis

## Data
    ORIG=$PROJ/original_data

## Get info from design file. Need a design file with unaligned reads from step 1 (aligning to junctions)

 DESIGN_FILE=$PROJ/design_files/design_file_alignment.csv
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}
    
## Set Different References
    REF=$PROJ/references/mm10_junctions_BT1

## Create Necessary Folders
    OUTPUT=$PROJ/aln_junctions
    if [ ! -e $OUTPUT ] ; then mkdir -p $OUTPUT; fi

    # Create JOB LOG directory and start log
        LOGS=$OUTPUT/job_logs 
        ALN_ERROR_LOG=$LOGS/size_errors.txt
        if [ ! -d $LOGS ]; then mkdir -p $LOGS; fi

        MYLOG=$LOGS/${NAME}.log
        printf "`date` $NAME PBS_ARRAYID:$PBS_ARRAYID HOSTNAME:$HOSTNAME\n" > $MYLOG
    # Create ALN LOG directory
        ALNLOGS=$OUTPUT/aln_logs 
        if [ ! -d $ALNLOGS ]; then mkdir -p $ALNLOGS; fi

    # Create UNALN READ directory
        UNALNDIR=$OUTPUT/unaln_reads
        if [ ! -d $UNALNDIR ]; then mkdir -p $UNALNDIR; fi

### FUNCTIONS FOR ALIGNMENT PIPELINE ###
# I have created a separate file that hold the alignment functions
source $PROJ/scripts/alignment_functions.sh

## Start Alignment Pipeline
    printf "<------------------- STARTING SE alignment process for $NAME [`date`] ------------------->\n" >> $MYLOG
    READS=$ORIG/${NAME}_SpliceMap.fq


    qual=`python /ufrc/mcintyre/share/python.git/identify_quality.py -i $READS`
    if [ $qual == 'phred64' ]; 
    then 
        # set to old illumina quality scores phred64/solexa 1.3
        btqual='--phred64-quals'
        lastqual='3'
    else
        # change to sanger format which is what all new illumina data is
        btqual='--phred33-quals'
        lastqual='1'
    fi

    bowtie_se_uniq

## Combine all Sam files 
    echo "START Combine SAM files">>$MYLOG
    cat *.sam >$OUTPUT/${NAME}.sam 2>>$MYLOG
    echo "FINISH Combining SAM files [`date`]" >>$MYLOG

## Combine all Unaln FASTQ 
   echo "START Combine Unaln FQ">>$MYLOG
   cat *_unaln_bt.fq >$UNALNDIR/${NAME}_unaln.fq 2>>$MYLOG
   echo "FINISH Combining SAM files [`date`]" >>$MYLOG

echo "Script complete, [`date`]" >>$MYLOG


