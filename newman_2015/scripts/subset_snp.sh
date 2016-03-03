#!/bin/bash

# Make a subset of the Onengut and eQTL SNPs

PROJ=/home/jrbnewman/concannon/eqtl_analysis
$PROJ/software/plink --bfile $PROJ/original_data/eurfam_03202014 --make-founders --extract $PROJ/pipeline_output/subset_snp_list.txt --make-bed --out $PROJ/original_data/eurfam_subset

