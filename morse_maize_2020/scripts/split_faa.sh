#!/bin/bash

module purge


## Split classification file so that isoproscan can be run in parallel


### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
PACBIO=/ufrc/mcintyre/share/maize_ainsworth/sqanti_post_filter
OUTPUT=$PROJ/split_faa
    mkdir -p ${OUTPUT}

## Set species
SPECIES=maize

## Get .faa file from SQANTI QC output
FAA=${PACBIO}/sqanti_filtered_corrected.faa

## Get number of lines to use in split
NUM=$(wc -l ${FAA} | awk '{if(int($1/20) % 2 == 1){ \
    print int($1/20)+1;} else{print int($1/20)+2}}')

## Split .faa file
split -d -l ${NUM} ${FAA} ${OUTPUT}/${SPECIES}.

## Make design file from list of split_faa files
DESIGN=$PROJ/design_files
if [ ! -e ${DESIGN} ]; then
    mkdir -p ${DESIGN}
fi
ls -1 ${OUTPUT} > ${DESIGN}/df_split_faa.csv
