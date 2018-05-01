#!/bin/sh

### Extract sequences for mm10 Refseq splicing catalog

PROJ=$MCLAB/conesa_pacbio

BED=$PROJ/created_files/conesa_refseq_mm10_splicing_catalogue_80bp.bed
FASTAREF=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
FASTAOUT=$PROJ/created_files/conesa_refseq_mm10_splicing_catalogue_80bp.tsv

bedtools getfasta -fi $FASTAREF \
                  -bed $BED \
                  -fo $FASTAOUT \
                  -name \
                  -split \
                  -tab \
                  -s

