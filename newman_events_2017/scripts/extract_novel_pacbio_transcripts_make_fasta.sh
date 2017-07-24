#!/bin/bash

## Subset only the novel PacBio transcripts for comparing against unannotated junctions

## Set directories

PROJ=$MCLAB/event_analysis

FASTAIN=$PROJ/references/Classification_of_CORRECTED_transcripts_Jeremy.txt
FASTAOUT=$PROJ/references/novel_pacbio_transcripts.fa

tr "\n" "\t" < $FASTAIN |
               sed -e 's/>/\n/g' |
               grep -v Known |
               sed -e 's/^/>/g' |
               tr "\t" "\n" | tail -n+2 > $FASTAOUT


