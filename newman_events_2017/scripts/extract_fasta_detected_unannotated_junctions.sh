#!/bin/bash

# Extract FASTA sequences of detected unannotated junctions

PROJ=$MCLAB/event_analysis
OUTPUT=$PROJ/references

INSEQ=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
BED=$OUTPUT/unannotated_junctions_detected_npc.bed
OUTSEQ=$OUTPUT/unannotated_junctions_detected_npc.fa

bedtools getfasta -fi $INSEQ -bed $BED -fo $OUTSEQ -split -name -s

