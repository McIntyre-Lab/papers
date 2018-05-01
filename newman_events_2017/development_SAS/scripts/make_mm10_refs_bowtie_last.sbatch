#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --partition=hpg1-compute
#SBATCH --job-name=make_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.make_refs.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb
#

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load bowtie/0.12.9
module load last/247

REFPATH=/ufrc/mcintyre/share/conesa_isoform_check/references
REF=mm10_for_bedtools_v2.fa
BTREF=mm10_for_bedtools_v2_BT1

cd $REFPATH

# make bowtie ref
bowtie-build $REF $BTREF

# make last ref
#lastdb $LASTREF $REF

