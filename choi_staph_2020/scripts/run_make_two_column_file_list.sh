#!/bin/bash


### create fq file list for input into stringMLST


PROJ=/home/ammorse/TB14/staph_relapse
INPUT=${PROJ}/SHARMA_6124_200205B6
OUTPUT=${PROJ}/stringMLST_analysis

python2 ${PROJ}/scripts/make_two_column_file_list_02dbs.py -d ${INPUT} -o ${OUTPUT}/read_file_list.txt -v


