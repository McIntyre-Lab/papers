#!/bin/bash

# Extract FASTA sequences of detected splicing

PROJ=$MCLAB/event_analysis
OUTPUT=$PROJ/references

INSEQ=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
BED=$OUTPUT/mm10_refseq_exon_fragments.bed
OUTSEQ=$OUTPUT/mm10_refseq_exon_fragments.fa

bedtools getfasta -fi $INSEQ -bed $BED -fo $OUTSEQ -name


