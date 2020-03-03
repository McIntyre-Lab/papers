#!/bin/bash

module purge


## Split classification file so that isoproscan can be run in parallel


### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
PACBIO=/ufrc/mcintyre/share/maize_ainsworth/sqanti_post_filter
OUTPUT=$PROJ/split_fasta
    mkdir -p ${OUTPUT}

## Set species
SPECIES=maize

## Get .fasta file from SQANTI QC output
## This is a multiline fasta file
FA=${PACBIO}/sqanti_filtered_corrected.fasta

## Convert to single line fasta
SL=${OUTPUT}/sqanti_filtered_corrected.single.fasta
awk '/^>/{printf("\n%s\n",$0);next;}{printf("%s",$0);} \
    END{printf("\n");}' ${FA} | tail -n +2 > ${SL}

## Get number of lines to use in split
NUM=$(wc -l ${SL} | awk '{if(int($1/20) % 2 == 1){ \
    print int($1/20)+1;} else{print int($1/20)+2}}')

## Split .fasta file
split -d -l ${NUM} ${SL} ${OUTPUT}/${SPECIES}.

## Remove single line fasta file
rm ${SL}

## Make design file from list of split_fasta files
DESIGN=$PROJ/design_files
if [ ! -e ${DESIGN} ]; then
    mkdir -p ${DESIGN}
fi
ls -1 ${OUTPUT} > ${DESIGN}/df_split_fasta.csv
