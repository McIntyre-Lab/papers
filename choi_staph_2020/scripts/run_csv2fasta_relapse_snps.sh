#!/bin/bash

## convert 2 column csv to fasta format
##  1st column contains header info and 2nd column contains seq

PROJ=/home/ammorse/mclab/SHARE/McIntyre_Lab/staph/seoung_ho_project

for REF in MSSA476 CA_347 ED98 Newman TCH60 ST20130941 
do

    sed 's/,/\n/g' $PROJ/output/fasta_snps_${REF}.csv > $PROJ/output/fasta_snps_${REF}.fasta

done

