#!/bin/bash

PROJ=/scratch/lfs/mcintyre/arbeitman_dsx-fru
HP=$PROJ/qc/splitdups_counts
OUT=$PROJ/qc/splitdups_summary.csv

cd $HP
FLAG=0

for FILE in *.csv
do 
	if [ $FLAG == 0 ] 
	then 
		cat $FILE > $OUT
		FLAG=1
	else
		tail -n +2 $FILE >> $OUT
	fi
done
