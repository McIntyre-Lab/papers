#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --partition=hpg1-compute
#SBATCH --account=mcintyre
#SBATCH --job-name=tophat
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 6:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --array=2-4

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export ROZ=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#export TMPDIR=/tmp

module load python
module load tophat/2.0.11
module load bowtie/0.12.9

## Set directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check

#    OUTPUT=$PROJ/tophat_junctions
#    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi



## Pull info from design file

    DESIGN_FILE=$PROJ/design_files/pacbio_samples.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
 
     READS=$PROJ/original_data/${SAMPLE}_SpliceMap.fq


     OUTPUT=$PROJ/tophat_junctions_${SAMPLE}
     if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi
 
BOWTIE_INDEXES=$PROJ/references

## Tophat requires Bowtie index files. Tophat will build them automatically if they do not exist. I will build them in the references folder.

    BTINDEX=${BOWTIE_INDEXES}/mm10_for_bedtools_v2_BT1

## Run TopHat

    tophat -o $OUTPUT --bowtie1 $BTINDEX $READS 

