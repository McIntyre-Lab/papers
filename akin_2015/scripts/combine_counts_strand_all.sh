#!/bin/bash
# Iterates through the forward and reverse strand coverage counts then combines them into one CSV for analysis

# Set directories

PROJ=/scratch/lfs/sugrue
FWD=$PROJ/coverage_fusions_sd
#REV=$PROJ/coverage_count_genes_rev

# Set output files

OUT=$PROJ/coverage_counts_fusions_strand.csv

# Combine FWD coverage

cd $FWD
FLAG=0
for FILE in *.csv;
     do
     if [ $FLAG == 0 ];
        then less $FILE > $OUT;
        FLAG=1;
     else
        tail -n+2 $FILE >> $OUT;
     fi;
done

# Combine REV coverage

#cd $REV
#for FILE in *.csv;
#     do
#     NAME=$(basename $FILE _rev.csv)
#     tail -n+2 $FILE | awk -v subj="$NAME," 'BEGIN{OFS=","} {print subj $0}' >> $OUT;
#     fi;
#done



