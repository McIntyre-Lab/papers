#!/bin/bash

## Concatenate mapped read counts for each sample mapped to each genome
## (B73, Mo17 CAU, and Mo17 YAN)
## Each genome is a separate file

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/check_isoseq3_lost_genes/minimap2

## Loop over genomes
for REF in b73 mo17_cau mo17_yan; do
    awk 'FNR==1 && NR!=1{next;}{print}' \
        ${IND}/*_2_${REF}_ccs_mapped_read_count.csv \
        > ${IND}/all_2_${REF}_ccs_mapped_read_count.csv
done
