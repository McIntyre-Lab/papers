#!/bin/bash

# Extract all transcript junctions from Ana's Pacbio transcripts
# to see what is common with what we can detect with the splicing
# catalog for mouse mm10 RefSeq

# Set folders

PROJ=$MCLAB/event_analysis
OUTPUT=$PROJ/analysis_output

# Location of script
SCRIPT=$MCLAB/junction_annotations/scripts/extractTranscriptJunctions.py

# Location of GFF3
GFF=$MCLAB/useful_mouse_data/mm10/gff/pacbio/conesa_pacbio_mm10.sorted_v2.gff

# Name of output
OUTNAME=mm10_conesa_pacbio_junctions.csv

# Run!
python $SCRIPT --gff $GFF --output $OUTPUT/$OUTNAME

