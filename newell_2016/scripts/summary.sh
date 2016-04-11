#!/bin/bash

PROJ=/home/fnew/mclab/arbeitman/arbeitman_ribotag
HP=$PROJ/pipeline_output/coverage_count_uniq_2
OUT=$HP/all_cvg_cnts_uniq.csv

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
