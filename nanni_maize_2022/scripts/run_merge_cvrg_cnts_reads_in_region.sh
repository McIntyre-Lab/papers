#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
DESIGN=$PROJ/design_files
INM=$PROJ/cvr_cnts_gene_mo17_cau
INB=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/cvr_cnts_shortRead_b73_gene
OUTM=${INM}
OUTB=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file

## Merge all reps in each genotype by gene
## Calculate sum of mapped reads by genotype_trt
## Output file for each genotype

## B73 mapped
python $PROJ/scripts/create_wide_cvrg_cnts_reads_in_region.py \
    -d ${INB} \
    -n primary_FBgn \
    -g B73,C123,Hp301,Mo17,NC338 \
    -p ${OUTB}/cvrg_shrtRead

## Mo17 CAU mapped
python $PROJ/scripts/create_wide_cvrg_cnts_reads_in_region.py \
    -d ${INM} \
    -n gene_id \
    -g B73,C123,Hp301,Mo17,NC338 \
    -p ${OUTM}/cvrg_shrtRead

