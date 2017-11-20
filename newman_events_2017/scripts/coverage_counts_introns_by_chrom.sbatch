#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=fusion_counts
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24g
#SBATCH --array=2-4
#
mkdir -p tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMPDIR=$(pwd)/tmp/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

module load python/2.7.6

    PROJ=/ufrc/mcintyre/share/conesa_isoform_check

#### Make Output Directory
        OUTPUT=$PROJ/coverage_counts_introns

echo "Making output directory"
        if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

## Design file
echo "Setting design file"
     DESIGN_FILE=$PROJ/design_files/pacbio_samples.csv

     DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
     IFS=',' read -ra ARRAY <<< "$DESIGN"

     SAMPLE=${ARRAY[0]}

     NAME=${SAMPLE}

## List of chromosomes including X and Y.
CHROM=/ufrc/mcintyre/share/conesa_isoform_check/references/chrom_list.txt


# Create LOG directory and start log

echo "Setting log"

        LOGS=$OUTPUT/logs
        if [ ! -e $LOGS ]; then mkdir -p $LOGS; fi
        
        MYLOG=$LOGS/${NAME}.log
        printf "`date` $FILE SLURM_ID:$SLURM_ARRAY_JOB_ID HOSTNAME:$HOSTNAME \n\n" > "${MYLOG}"
echo "Log set as $MYLOG "

#### COVERAGE COUNTS
    #BED=$PROJ/references/conesa_pacbio_mm10_fusions_si_coverage.bed
    BED=$PROJ/references/mm10_introns_from_fusions_calculated.bed
    SAM=$PROJ/bwa_mem_aln_pe/bwa_split/${NAME}_mapped_all.sam
    MPILEUP=$PROJ/genome_mpileups/${NAME}.mpileup


echo "Coverage counts"

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
