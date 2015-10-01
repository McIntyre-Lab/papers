#!/bin/bash
#PBS -M jfear@ufl.edu
#PBS -N sampy
#PBS -r n
#PBS -q bio
#PBS -o /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/ase_pipeline/
#PBS -e /scratch/lfs/mcintyre/cegs_ase_paper/scripts/PBS_LOGS/ase_pipeline/
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=35gb
##PBS -t 1-417
#PBS -t 51

set -o nounset
set -e

module load python/2.7.6

# Set Directories
PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
INDIR=$PROJ/fastq_split
REF=/scratch/lfs/mcintyre/references/dmel_fb551/dmel-all-chromosome-r5.51.fasta
FUS=/scratch/lfs/mcintyre/references/dmel_fb551/dmel-non-redundant-r5.51_fusions.bed

# Set Design info
DESIGN_FILE=$PROJ/design_files/CEGS_68_lines_no_tech.txt
DESIGN=$(sed -n "${PBS_ARRAYID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

LINE=${ARRAY[0]}
MV=${ARRAY[1]}
REP=${ARRAY[2]}
NAME=${LINE}_${MV}${REP}

# Create output directories and logs
OUTDIR=$PROJ/ase_pipeline_output/ase_counts_fb551_updated_fusions
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

LOGS=$OUTDIR/logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi

MYLOG=$LOGS/${NAME}.log
printf "`date` $LINE $MV $REP SGE_ID:$PBS_ARRAYID HOSTNAME:$HOSTNAME \n\n" > $MYLOG

# Combine reads for sam-compare
READS=$TMPDIR/${NAME}_distinct.fq
cat $INDIR/${NAME}.*_distinct.fq > $READS 

# Grab Sam files
SAM1=$PROJ/ase_pipeline_output/aln_upd_fusions/${NAME}_${LINE}_upd_fusions.sam
SAM2=$PROJ/ase_pipeline_output/aln_upd_fusions/${NAME}_w1118_upd_fusions.sam

python /scratch/lfs/jfear/devel/sam-compare/sam_compare.py \
    -l 95 \
    -f $FUS \
    -q $READS \
    -A $SAM1 \
    -B $SAM2 \
    -c $OUTDIR/ase_counts_${NAME}.csv \
    -t $OUTDIR/ase_totals_${NAME}.txt \
    -g $LOGS/ase_counts_${NAME}.log 

echo "Script Complete, [`date`]" >>$MYLOG
