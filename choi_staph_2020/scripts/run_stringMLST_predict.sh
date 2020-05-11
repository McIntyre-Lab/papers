#!/bin/bash 

## stringMLST Predict

## Default kmer of 35


PROJ=/home/ammorse/TB14/staph_relapse

cd $PROJ/stringMLST_analysis


stringMLST.py --predict \
    -l read_file_list.txt \
    -P staphDB/Staphylococcus_aureus \
    -o stringMLST_out.tsv \
    -x

