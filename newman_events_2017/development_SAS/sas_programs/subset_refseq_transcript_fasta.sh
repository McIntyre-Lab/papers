#!/bin/sh

### Subset RefSeq transcript FASTA for RSEM

PROJ=$MCLAB/event_analysis

LIST1=refseq_list_exp_transcripts.txt
LIST2=refseq_list_exp_transcripts_75perc_dtct.txt
LIST3=refseq_list_exp_transcripts_100perc_dtct.txt

SCRIPT=$PROJ/scripts/subset_fasta_file_by_ID_list.py

FASTAIN=$PROJ/references/refseq_mm10.fa
FASTAOUT1=$PROJ/references/refseq_mm10_exp_transcripts.fa
FASTAOUT2=$PROJ/references/refseq_mm10_exp_transcripts_75perc_dtct.fa
FASTAOUT3=$PROJ/references/refseq_mm10_exp_transcripts_100perc_dtct.fa


python2.7 $SCRIPT -i $FASTAIN -l $LIST1 -o $FASTAOUT1
python2.7 $SCRIPT -i $FASTAIN -l $LIST2 -o $FASTAOUT2
python2.7 $SCRIPT -i $FASTAIN -l $LIST3 -o $FASTAOUT3

