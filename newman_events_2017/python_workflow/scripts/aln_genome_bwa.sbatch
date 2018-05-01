#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=aln_genome_bwa
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/event_analysis/scripts/SLURM_LOGS/out.aln_bwa.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16gb
#SBATCH --array=1-4
#

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load bwa/0.7.12
module load python/2.7.6
module load samtools

## Set directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check/event_analysis
    ORIG=$PROJ/original_data
    
    OUTPUT=$PROJ/bwa_mem_aln_pe
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

    ## Output folder for split SAM files
    SPLITSAM=$OUTPUT/bwa_split
    if [ ! -e $SPLITSAM ]; then mkdir -p $SPLITSAM; fi

    ## Design file
    
    
 DESIGN_FILE=$PROJ/design_files/design_file_alignment.csv
 DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
 IFS=',' read -ra ARRAY <<< "$DESIGN"

 SAMPLE=${ARRAY[0]}
 NAME=${SAMPLE}

## set references
    REF=$PROJ/references/mm10_for_bedtools_v2_BWA

#### create bwa-mem reference

#    bwa index -p $REF -a bwtsw $FASTA

#### run bwa-mem

    bwa mem -t 4 -M $REF $ORIG/${NAME}_SpliceMap.fq > $OUTPUT/${NAME}.sam

#### parse bwa alignment file
    python $PROJ/scripts/BWASplitSAM_07mai.py \
           -fq1 $ORIG/${NAME}_SpliceMap.fq \
            -s $OUTPUT/${NAME}.sam \
            --outdir $SPLITSAM

