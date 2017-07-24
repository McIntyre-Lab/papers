#!/bin/bash

## Make BLAST DB for aligning unannotated junctions and IR events against Pacbio transcripts


# Set directories
PROJ=$MCLAB/event_analysis
FASTA=$PROJ/references/Classification_of_CORRECTED_transcripts_Jeremy.txt
OUTPUT=$PROJ/references

if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT ; fi

OUTNAME=pacbio_transcripts.BLAST

# Make BLAST DB

makeblastdb -in $FASTA -out $OUTPUT/$OUTNAME -dbtype nucl
