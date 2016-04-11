#!/bin/bash

MAST=/home/jfear/opt/meme/bin/mast
PROJ=$MCLAB/arbeitman_fru_network
REF=$MCLAB/useful_dmel_data/flybase530/dmel_sequence/dmel-all-chromosome-r5.30.fasta

for motif in a b c
do
	$MAST $PROJ/motif_analysis/meme_fru_${motif}/meme.xml $REF -o $PROJ/motif_analysis/mast_fru_${motif} &>$PROJ/motif_analysis/mast_fru_${motif}.log 
done
