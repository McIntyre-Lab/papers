#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --partition=hpg1-compute
#SBATCH --account=mcintyre
#SBATCH --job-name=rsem_calc
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16gb
#SBATCH --array=17-21
#

## 20




mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export ROZ=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}



module load R
module load bowtie/0.12.9
module load rsem/1.2.28

PROJ=/ufrc/mcintyre/share/conesa_isoform_check
REF=$PROJ/references/rsem
OUTDIR=$PROJ/rsem_output

if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi



    DESIGN_FILE=$PROJ/design_files/rsem_pacbio_design2.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    REFSET=${ARRAY[1]}
   
    RSEMREF=${REFSET}_RSEM
    
    OUTPUT=$OUTDIR/${REFSET}_${SAMPLE}
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT ; fi


     READS=$PROJ/original_data/${SAMPLE}_SpliceMap.fq



## Run RSEM

rsem-calculate-expression -p 4 \
                          --temporary-folder ${ROZ} \
                          --calc-ci \
                          --estimate-rspd \
                          --fragment-length-mean 80 \
                          --fragment-length-sd 50 \
                          ${READS} ${REF}/${RSEMREF} ${OUTPUT}

