#!/bin/bash


## Folder of split gff3 files

INGFFS=../gff/gff3/split_gffs
OUTGFF=../gff/gff3/aceview.gff

FLAG=0
for FILE in $INGFFS;
do
   if [ $FLAG == 0 ];
   then less $FILE > $OUTGFF;
   FLAG=1;
   else
   tail -n +2 $FILE >> $OUTGFF;
   fi;
done
   
