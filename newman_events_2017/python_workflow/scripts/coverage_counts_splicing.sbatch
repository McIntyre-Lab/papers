#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=splicing_counts
#SBATCH --mail-type=FAIL
#SBATCH --partition=hpg1-compute
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/event_analysis/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --array=1-4
#
mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load python/2.7.6

    PROJ=/ufrc/mcintyre/share/conesa_isoform_check/event_analysis

#### Make Output Directory
        OUTPUT=$PROJ/coverage_counts_splicing

        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

## Design file
     DESIGN_FILE=$PROJ/design_files/design_file_alignment.csv

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}

     NAME=${SAMPLE}


# Create LOG directory and start log
        LOGS=$OUTPUT/logs
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        MYLOG=$LOGS/${FILE}.log
        printf "`date` $FILE SLURM_ID:$SLURM_ARRAY_JOB_ID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}"

#### COVERAGE COUNTS
    #BED=$PROJ/references/conesa_pacbio_mm10_splicing_catalogue_80bp_for_coverage.bed
    BED=$PROJ/references/mm10_junctions_coverage.bed
    SAM=$PROJ/aln_junctions/${NAME}.sam
    MPILEUP=$PROJ/splicing_mpileups/${NAME}.mpileup

    echo "Starting Coverage Counts for $NAME (plus strand) `date`" > "${MYLOG}"
    python $PROJ/scripts/rpkm_calculate.py \
        -b $BED \
        -m $MPILEUP \
        -s $SAM \
        -n ${NAME} \
        --cv \
        -g "${MYLOG}" \
        -o $OUTPUT/cvrg_cnts_${NAME}.csv
    echo "Finished Coverage Counts for $NAME `date`" >> "${MYLOG}"

echo "Script complete [`date`]" >> "${MYLOG}"

