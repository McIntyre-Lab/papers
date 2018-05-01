#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=bwa_index
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/event_analysis/scripts/SLURM_LOGS/out.bwa_index.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#


mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}


module load bwa/0.7.12 samtools/1.3.1

## Set Directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check/event_analysis/references
    REF=mm10_for_bedtools_v2.fa

cd $PROJ
samtools faidx $REF
bwa index -a bwtsw $PROJ/$REF -p mm10_for_bedtools_v2_BWA

