#!/bin/bash

## BLAST detected splicing events against PacBio transcript sequences and RefSeq transcript sequences
## I want to double-check that events are BLASTing where they should!!!

PROJ=$MCLAB/event_analysis

INPUT=$PROJ/references/refseq_mm10_detected_splicing_events_nsc.fa

PACBIODB=$PROJ/references/pacbio_transcripts.BLAST
REFSEQDB=$PROJ/references/refseq_transcripts.BLAST

OUTDIR=$PROJ/analysis_output/blast_output

if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

OUTPUT1=$OUTDIR/blast_events_to_pacbio_basic.tsv
OUTPUT2=$OUTDIR/blast_events_to_pacbio_min_length_50.tsv

OUTPUT3=$OUTDIR/blast_events_to_refseq_basic.tsv
OUTPUT4=$OUTDIR/blast_events_to_refseq_min_length_50.tsv

### PACBIO BLAST

## Basic BLAST
#blastn -db $PACBIODB -num_threads 1 -query $INPUT -outfmt 6 > $OUTPUT1

## Minimum match size is 50bp and 90% matching. This should filter out any match that only matches the donor exonic sequence
blastn -db $PACBIODB -num_threads 1 -query $INPUT -outfmt 6 -task megablast -dust no -strand plus -evalue 1 -word_size 50 -perc_identity 90 > $OUTPUT2


### REFSEQ BLAST

## Basic BLAST
#blastn -db $REFSEQDB -num_threads 1 -query $INPUT -outfmt 6 > $OUTPUT3

## Minimum match size is 50bp and 90% matching. This should filter out any match that only matches the donor exonic sequence
blastn -db $REFSEQDB -num_threads 1 -query $INPUT -outfmt 6 -task megablast -dust no -strand plus -evalue 1 -word_size 50 -perc_identity 90 > $OUTPUT4

