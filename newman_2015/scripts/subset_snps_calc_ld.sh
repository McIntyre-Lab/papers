#!/bin/bash

# Calculate LD between Onengut and eQTL SNPs

PROJ=/home/jrbnewman/concannon/eqtl_analysis
$PROJ/software/plink2 --bfile $PROJ/original_data/eurfam_subset --r2 --ld-window-kb 10000 --ld-window-r2 0 --ld-snp-list $PROJ/pipeline_output/eqtl_snp_list_1.txt --out $PROJ/pipeline_output/eurfam_subset_ld_1

$PROJ/software/plink2 --bfile $PROJ/original_data/eurfam_subset --r2 --ld-window-kb 10000 --ld-window-r2 0 --ld-snp-list $PROJ/pipeline_output/eqtl_snp_list_2.txt --out $PROJ/pipeline_output/eurfam_subset_ld_2

## Process putput so it's easier to import into SAS

awk '{print $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7}' $PROJ/pipeline_output/eurfam_subset_ld_1.ld > $PROJ/pipeline_output/eurfam_subset_ld.csv
tail -n+2 $PROJ/pipeline_output/eurfam_subset_ld_2.ld | awk '{print $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7}' >> $PROJ/pipeline_output/eurfam_subset_ld.csv


