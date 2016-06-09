#!/bin/bash
#PBS -M jfear@ufl.edu
#PBS -N alnmsk310369
#PBS -r n
#PBS -q bio
#PBS -o /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/aln_2_masked/
#PBS -e /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/aln_2_masked/
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=2
#PBS -l pmem=15gb
#PBS -t 1-417

module load bowtie/0.12.9 
module load last/247 
module load python/2.7.6
module load samtools/1.1

# Store number of processors used -- i.e. number of files to split into to run LAST
NUMPROCS=2

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
INDIR=$PROJ/fastq_split

# Set Design info
DESIGN_FILE=$PROJ/design_files/CEGS_68_lines_no_tech.txt
DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

LINE=${ARRAY[0]}
MV=${ARRAY[1]}
REP=${ARRAY[2]}
NAME=${LINE}_${MV}${REP}

# Set Different References
REF=$PROJ/references/masked_genome/${LINE}/w1118_2_${LINE}_masked_genome_BT1
LASTREF=$PROJ/references/masked_genome/${LINE}/w1118_2_${LINE}_masked_genome_LAST
FAREF=$PROJ/references/masked_genome/${LINE}/w1118_2_${LINE}_masked_genome.fasta

# Create output directories and logs
OUTDIR=$PROJ/ase_pipeline_output/aln_masked_genome
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

LOGS=$OUTDIR/logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi

MYLOG=$LOGS/${NAME}_masked.log
printf "`date` $LINE $MV $REP SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

ALNLOGS=$OUTDIR/aln_logs 
if [ ! -e $ALNLOGS ]; then mkdir -p $ALNLOGS; fi

# Start Alignment Pipeline
## FUNCTIONS FOR ALIGNMENT PIPELINE ###
source $PROJ/scripts/alignment_functions.sh

## Start Alignment Pipeline
    printf "<------------------- STARTING SE alignment process for $NAME [`date`] ------------------->\n" >> $MYLOG
    READS=$TMPDIR/${NAME}_distinct.fq
    cat $INDIR/${LINE}_${MV}${REP}.*_distinct.fq > $READS 

    qual=`python /scratch/lfs/mcintyre/python.git/identify_quality.py -i $READS`
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

    bowtie_se_all
    last_se_all

## Combine all Sam files
    echo "START Combine SAM files">>$MYLOG
    cat *.sam >$TMPDIR/${NAME}.sam 2>>$MYLOG
    echo "FINISH Combining SAM files [`date`]" >>$MYLOG

## Convert to BAM files
    echo "Convert SAM to BAM">>$MYLOG
    samtools view -bt $FAREF.fai -o $TMPDIR/$NAME.bam $TMPDIR/$NAME.sam 2>>$MYLOG
    samtools sort -m 10000000000 $TMPDIR/$NAME.bam $TMPDIR/$NAME.sorted 2>>$MYLOG
    samtools index $TMPDIR/$NAME.sorted.bam 2>>$MYLOG
    echo "Finished converting SAM to BAM">>$MYLOG
    cd $TMPDIR
    cp $NAME.sorted.* $OUTDIR/

echo "Script Complete, [`date`]" >>$MYLOG
