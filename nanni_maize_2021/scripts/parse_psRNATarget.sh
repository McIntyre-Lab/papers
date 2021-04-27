#!/bin/bash

module purge


## Parse psRNATarget output file

### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
MIRNA=$PROJ/miRNA_target

## Set input and output files
PSRNA=${MIRNA}/psRNATargetJob-*.txt
PARSE=${MIRNA}/parsed_psRNATarget.txt

awk '{if(NR>2){print $2"\t"$7"\t"$8"\tID="$1";Note="$11}}' \
    ${PSRNA} > ${PARSE}
