#!/bin/sh

## Cat all parsed same file summaries from bwa pe and se separately (have different headers)

PROJ=/ufrc/mcintyre/share/etoh_srna
IN=$PROJ/rnaseq/aln_genome/bwa_parsed_sam_files

awk 'FNR==1 && NR!=1{next;}{print}' $IN/*_*_bwa_SE_summary.csv \
    > $IN/combined_bwa_SE_summary.csv


awk 'FNR==1 && NR!=1{next;}{print}' $IN/*_*_bwa_PE_summary.csv \
    > $IN/combined_bwa_PE_summary.csv
