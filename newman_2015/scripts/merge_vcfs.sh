#!/bin/bash

FLAG=0
for FILE in ../immunobase_snps/*.vcf
   do
   if [ $FLAG == 0 ];
      then less $FILE > ../pipeline_output/merged_vcf.tsv;
      FLAG=1;
   else
      less $FILE >> ../pipeline_output/merged_vcf.tsv;
   fi
done

