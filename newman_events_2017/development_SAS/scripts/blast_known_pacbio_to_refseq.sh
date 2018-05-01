#!/bin/bash

#### BLAST known PacBio transcripts against RefSeq transcripts
### Using this to identify the corresponding RefSeq IDs for each PacBio transcript

PROJ=$MCLAB/event_analysis
REF=$PROJ/references

QUERY=$REF/known_pacbio_transcripts.fa
DB=$REF/refseq_transcripts.BLAST

OUT=$PROJ/analysis_output/blast_output/known_pacbio2refseq_best_hits.tsv

## Run MegaBLAST

blastn -db $DB \
       -num_threads 1 \
       -query $QUERY \
       -outfmt 6 \
       -task megablast \
       -evalue 1e-10 \
       -best_hit_score_edge 0.05 \
       -best_hit_overhang 0.25 \
       -perc_identity 80 \
       -max_target_seqs 1 > $OUT

