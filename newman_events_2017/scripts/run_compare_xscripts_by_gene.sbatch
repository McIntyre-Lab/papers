#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=xs_pairs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --array=2-30
#

##30

mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load sas/9.4

#Set up variables

PROJ=/ufrc/mcintyre/share/conesa_isoform_check
SCRIPTS=$PROJ/sas_programs
LOGS=$PROJ/sas_programs/sas_logs
if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
SASDATA=$PROJ/sas_data
if [ ! -e $SASDATA ]; then mkdir -p $SASDATA; fi


# Pull info from design file
    DESIGN_FILE=$PROJ/design_files/pbiso_xscript_compare_design_file.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    STARTGENE=${ARRAY[0]}
    STOPGENE=${ARRAY[1]}
    SETNUM=${ARRAY[2]}

## Set SAS program to run
PRG0=analysis_complete_isoforms_compare_pairwise_jrbn2.sas

## set other variables:
GROUP_PREFIX=pbiso
SASLIB=$PROJ/sas_data


SASTEMP=/ufrc/mcintyre/share/conesa_isoform_check/sas_temp/sas_work_${SLURM_ARRAY_TASK_ID}
if [ ! -e $SASTEMP ]; then mkdir -p $SASTEMP; fi
cd $SASDATA

sas -log $LOGS/${PRG0}_${SLURM_ARRAY_TASK_ID}.log -work $SASTEMP -sysin $SCRIPTS/$PRG0 -set geneStart ${STARTGENE} -set geneStop ${STOPGENE} -set setNum ${SETNUM} -set libpath $SASLIB -set group ${GROUP_PREFIX} -memsize 20G 

rm $TMPDIR/*

