#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=intron_counts
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/event_analysis/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36gb
#SBATCH --array=1-4
#
mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load python/2.7.6

    PROJ=/ufrc/mcintyre/share/conesa_isoform_check/event_analysis

#### Make Output Directory
        OUTPUT=$PROJ/coverage_counts_introns

        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

## Design file
     DESIGN_FILE=$PROJ/design_files/design_file_alignment.csv

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}

     NAME=${SAMPLE}

## List of chromosomes including X and Y.
CHROM=/ufrc/mcintyre/share/conesa_isoform_check/references/chrom_list.txt

# Create LOG directory and start log
        LOGS=$OUTPUT/logs
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        MYLOG=$LOGS/${FILE}.log
        printf "`date` $FILE SLURM_ID:$SLURM_ARRAY_JOB_ID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}"

#### COVERAGE COUNTS
    BED=$PROJ/references/mm10_introns_from_fusions.bed
    SAM=$PROJ/bwa_mem_aln_pe/${NAME}.sam
    MPILEUP=$PROJ/genome_mpileups/${NAME}.mpileup


   echo "Starting Coverage Counts for $NAME (plus strand) `date`" > "${MYLOG}"

COUNT=0
for chrom in $(cat $CHROM)
do

awk -v chr="$chrom" '{if ($3 == chr) print $0}'

    awk -v chr="$chrom" '{if ($3 == chr) print $0}' $SAM > $TMPDIR/${NAME}_${chrom}.sam
    awk -v chr="$chrom" '{if ($1 == chr) print $0}' $BED > $TMPDIR/mm10_introns_${chrom}.bed

    echo "Starting Coverage Counts for $NAME (plus strand) `date`" > "${MYLOG}"
    python $PROJ/scripts/rpkm_calculate.py \
        -b $TMPDIR/mm10_introns_${chrom}.bed \
        -m $MPILEUP \
        -s $TMPDIR/${NAME}_${chrom}.sam \
        -n ${NAME} \
        --cv \
        -g "${MYLOG}" \
        -o $TMPDIR/${NAME}_${chrom}.csv
    echo "Finished Coverage Counts for $NAME `date`" >> "${MYLOG}"

   if [ $COUNT == 0 ]
   then
      cat $TMPDIR/${NAME}_${chrom}.csv >  $OUTPUT/${NAME}.csv
        COUNT=1
    else
        tail -n+2 $TMPDIR/${NAME}_${chrom}.csv >> $OUTPUT/${NAME}.csv
    fi
done

echo "Script complete [`date`]" >> "${MYLOG}"

echo "done"
