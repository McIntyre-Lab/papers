#!/bin/bash

## BLAST detected splicing events against PacBio transcript sequences

PROJ=$MCLAB/event_analysis

INPUT=$PROJ/references/refseq_mm10_detected_splicing_events_nsc.fa

PACBIODB=$PROJ/references/pacbio_transcripts.BLAST

OUTDIR=$PROJ/analysis_output/blast_output

if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

OUTPUT1=$OUTDIR/blast_events_to_pacbio_basici_nomulti.tsv
OUTPUT2=$OUTDIR/blast_events_to_pacbio_min_length_50_nomulti.tsv


### PACBIO BLAST

## Basic BLAST
blastn -db $PACBIODB -num_threads 1 -query $INPUT -outfmt 6 > $OUTPUT1

## Minimum match size is 50bp and 90% matching. This should filter out any match that only matches the donor exonic sequence
blastn -db $PACBIODB -num_threads 1 -query $INPUT -outfmt 6 -task megablast -dust no -strand plus -evalue 1 -word_size 50 -perc_identity 90 > $OUTPUT2


