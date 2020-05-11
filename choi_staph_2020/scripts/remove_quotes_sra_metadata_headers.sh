#!/bin/bash 

## remove quotes from sra design file metadata table headers

PROJ=/home/ammorse/staph_relapse/sra

sed -i 's/"//g' $PROJ/design_sra_meta.tsv


