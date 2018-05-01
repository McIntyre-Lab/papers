#!/bin/bash

## Combine fusion coverage counts

PROJ=/ufrc/mcintyre/share/conesa_isoform_check

FUSCOUNTS=$PROJ/coverage_counts_splicing
OUTCOUNTS=$PROJ/mm10_refseq_splicing_counts.csv

cd $FUSCOUNTS

FLAG=0
for FILE in *.csv;
do
   if [ $FLAG == 0 ]
   then cat $FILE > $OUTCOUNTS
   FLAG=1
   else
   tail -n+2 $FILE >> $OUTCOUNTS
   fi
done

