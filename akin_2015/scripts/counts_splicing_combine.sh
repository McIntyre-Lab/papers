#!/bin/bash

# Iterates through each sample's batch of CSVs and merges them together

MDIR=/scratch/lfs/counts_hg19_splicing/merged
mkdir $MDIR
num=1
while [ $num -lt 65 ];
   do
   NAME=$(awk -v "line=$num" 'NR == line' /scratch/lfs/sugrue/design_files/sample_list.txt)
   FLAG=0
   MERGE=$MDIR/${NAME}.csv
   for FILE in ${NAME}*;
        do
        if [ $FLAG == 0 ];
     	   then less $FILE > $MERGE;
           FLAG=1;
        else
           tail -n+2 $FILE >> $MERGE;
        fi;
   done
   num=$[$num+1]
done

# Takes merged CSVs, appends sample_id to each and outputs as one CSV

OUT=/scratch/lfs/sugrue/coverage_counts_splicing.csv
echo sample,fusion_id,mapped_reads,read_length,region_length,region_depth,reads_in_region,apn,rpkm,mean,std,cv > $OUT

cd $MDIR
for FILE in *.csv;
     do
     NAME=$(basename $FILE .csv)
     SUBJECT=$( echo $NAME | sed -e 's/cvrg_cnts_//g' )
     tail -n+2 $FILE | awk -v subj="$SUBJECT," 'BEGIN{OFS=","} {print subj $0}' >> $OUT;
done

