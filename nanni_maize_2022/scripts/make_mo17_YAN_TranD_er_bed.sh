#!/bin/bash

## Make BED file of Mo17 YAN fusions (exonic regions) from TranD output

## Set directories
REF=/blue/mcintyre/share/references/maize_mo17
IND=${REF}/TranD_gene_maize_mo17_YAN

## Columns:
## er_id    2
## er_chr   3
## er_start 4
## er_end   5
## star coord is 0-based like a bed file so no changing required for coords
awk -F "," 'NR!=1{print $3"\t"$4"\t"$5"\t"$2}' \
    ${IND}/event_analysis_er.csv \
    > ${REF}/mo17_YAN_TranD_exon_region_coverage.bed
