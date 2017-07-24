#!/bin/bash

# Extract FASTA sequences of reclassified IR events

PROJ=$MCLAB/event_analysis
OUTPUT=$PROJ/references

INSEQ=$MCLAB/useful_mouse_data/mm10/fasta/mm10_for_bedtools_v2.fa
BED=$OUTPUT/reclassified_ir_events_npc.bed
OUTSEQ=$OUTPUT/reclassified_ir_events_npc.fa


bedtools getfasta -fi $INSEQ -bed $BED -fo $OUTSEQ -split -name -s

