#!/bin/bash

## Extract RefSeq transcript sequences from genome FASTA

## Set directories

PROJ=$MCLAB/event_analysis

REF=$PROJ/references

### (1) Make BED12 file for extracting sequences

GFF=$MCLAB///gff/pacbio/aconesa_refseq_gff3_v4.gff

BEDOUT=$REF/RefSeq.bed

python $PROJ/scripts/make_transcript_bed.py --gff $GFF --output $BEDOUT


#### (2) Remove erroneous transcript ID
## There is an erroneous transcript consisting of zero exons in the GFF3 file
## I want to remove this so that bedtools doesn't throw back an error

grep -v NM_009081_A $BEDOUT > $REF/RefSeq2.bed


### (3) Extract FASTA sequences

FASTAIN=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
FASTAOUT=$REF/refseq_mm10.fa

bedtools getfasta -name -s -split -fi $FASTAIN -bed $REF/RefSeq2.bed -fo $FASTAOUT

