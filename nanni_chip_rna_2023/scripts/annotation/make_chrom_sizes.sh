#!/bin/bash

## Generate chrom.sizes file if it does not already exist

## Require global variable FASTA to be set to the FASTA of interest
## Generates FASTA index file (.fai) if it does not already exist
## Outputs .chrom.sizes file with the same prefix as the FASTA

CHROM=$(dirname ${FASTA})/$(basename ${FASTA} .fasta).chrom.sizes

if [ ! -e ${CHROM} ]; then
    if [ ! -e ${FASTA}.fai ]; then
        samtools faidx ${FASTA}
    fi
    cut -f 1,2 ${FASTA}.fai > ${CHROM}
fi
