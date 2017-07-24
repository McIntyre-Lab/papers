#!/bin/bash

# Extract FASTA sequences of detected splicing

PROJ=$MCLAB/event_analysis
OUTPUT=$PROJ/references

INSEQ=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
BED=$OUTPUT/refseq_mm10_detected_splicing_events_nsc_nomulti.bed
OUTSEQ=$OUTPUT/refseq_mm10_detected_splicing_events_nsc_nomulti.fa

bedtools getfasta -fi $INSEQ -bed $BED -fo $OUTSEQ -split -name -s


