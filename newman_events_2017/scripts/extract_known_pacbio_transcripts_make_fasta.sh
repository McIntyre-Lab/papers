#!/bin/bash

## Subset only the known PacBio transcripts for comparing against RefSeq transcripts

## Set directories

PROJ=$MCLAB/event_analysis

FASTAIN=$PROJ/references/Classification_of_CORRECTED_transcripts_Jeremy.txt
FASTAOUT=$PROJ/references/known_pacbio_transcripts.fa

tr "\n" "\t" < $FASTAIN |
               sed -e 's/>/\n/g' |
               grep Known |
               sed -e 's/^/>/g' |
               tr "\t" "\n" > $FASTAOUT


