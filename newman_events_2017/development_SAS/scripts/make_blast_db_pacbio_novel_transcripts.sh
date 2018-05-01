#!/bin/bash

## Make BLAST DB for aligning unannotated junctions and IR events against Pacbio transcripts


# Set directories
PROJ=$MCLAB/event_analysis
FASTA=$PROJ/references/novel_pacbio_transcripts.fa
OUTPUT=$PROJ/references

if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT ; fi

OUTNAME=novel_pacbio_transcripts.BLAST

# Make BLAST DB

makeblastdb -in $FASTA -out $OUTPUT/$OUTNAME -dbtype nucl
