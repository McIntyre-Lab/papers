#!/bin/bash

## BLAST reclassified IR events against novel PacBio transcript sequences


PROJ=$MCLAB/event_analysis

INPUT=$PROJ/references/reclassified_ir_events_npc.fa

DB=$PROJ/references/pacbio_transcripts.BLAST

OUTDIR=$PROJ/analysis_output/blast_output

if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

OUTPUT1=$OUTDIR/blast_reclassified_ir_events_to_pacbio_basic.tsv
OUTPUT2=$OUTDIR/blast_reclassified_ir_events_to_pacbio_min_length_50.tsv

## Minimum match size is 50bp and 90% matching. This should filter out any match that only matches the donor exonic sequence
blastn -db $DB -num_threads 1 -query $INPUT -outfmt 6 -task megablast -dust no -strand plus -evalue 1 -word_size 50 -perc_identity 90 > $OUTPUT1

## Basic BLAST
#blastn -db $DB -num_threads 1 -query $INPUT -outfmt 6 > $OUTPUT2

