#!/bin/sh

### Subset RefSeq transcript FASTA for RSEM

PROJ=$MCLAB/event_analysis

#LIST1=$PROJ/analysis_output/refseq_list_exp_transcripts.txt
#LIST2=$PROJ/analysis_output/refseq_list_exp_transcripts_75perc_dtct.txt
#LIST3=$PROJ/analysis_output/refseq_list_exp_transcripts_100perc_dtct.txt
#LIST4=$PROJ/analysis_output/refseq_list_exp_transcripts_75perc_dtct_apn5.txt
#LIST5=$PROJ/analysis_output/refseq_list_exp_transcripts_75perc_dtct_apn10.txt

LIST6=$PROJ/analysis_output/refseq_list_exp_transcripts_50perc_dtct.txt
LIST7=$PROJ/analysis_output/refseq_list_exp_transcripts_50perc_dtct_apn5.txt
LIST8=$PROJ/analysis_output/refseq_list_exp_transcripts_50perc_dtct_apn10.txt
LIST9=$PROJ/analysis_output/refseq_list_exp_transcripts_100perc_dtct_apn5.txt
LIST10=$PROJ/analysis_output/refseq_list_exp_transcripts_100perc_dtct_apn10.txt


SCRIPT=$PROJ/scripts/subset_fasta_file_by_ID_list.py

FASTAIN=$PROJ/references/refseq_mm10.fa
#FASTAOUT1=$PROJ/references/refseq_mm10_exp_transcripts.fa
#FASTAOUT2=$PROJ/references/refseq_mm10_exp_transcripts_75perc_dtct.fa
#FASTAOUT3=$PROJ/references/refseq_mm10_exp_transcripts_100perc_dtct.fa
#FASTAOUT4=$PROJ/references/refseq_mm10_exp_transcripts_75perc_dtct_apn5.fa
#FASTAOUT5=$PROJ/references/refseq_mm10_exp_transcripts_75perc_dtct_apn10.fa

FASTAOUT6=$PROJ/references/refseq_mm10_exp_transcripts_50perc_dtct.fa
FASTAOUT7=$PROJ/references/refseq_mm10_exp_transcripts_50perc_dtct_apn5.fa
FASTAOUT8=$PROJ/references/refseq_mm10_exp_transcripts_50perc_dtct_apn10.fa
FASTAOUT9=$PROJ/references/refseq_mm10_exp_transcripts_100perc_dtct_apn5.fa
FASTAOUT10=$PROJ/references/refseq_mm10_exp_transcripts_100perc_dtct_apn10.fa


#python2.7 $SCRIPT -i $FASTAIN -l $LIST1 -o $FASTAOUT1
#python2.7 $SCRIPT -i $FASTAIN -l $LIST2 -o $FASTAOUT2
#python2.7 $SCRIPT -i $FASTAIN -l $LIST3 -o $FASTAOUT3
#python2.7 $SCRIPT -i $FASTAIN -l $LIST4 -o $FASTAOUT4
#python2.7 $SCRIPT -i $FASTAIN -l $LIST5 -o $FASTAOUT5

python2.7 $SCRIPT -i $FASTAIN -l $LIST6 -o $FASTAOUT6
python2.7 $SCRIPT -i $FASTAIN -l $LIST7 -o $FASTAOUT7
python2.7 $SCRIPT -i $FASTAIN -l $LIST8 -o $FASTAOUT8
python2.7 $SCRIPT -i $FASTAIN -l $LIST9 -o $FASTAOUT9
python2.7 $SCRIPT -i $FASTAIN -l $LIST10 -o $FASTAOUT10

SCRIPT=$PROJ/scripts/subset_fasta_file_by_ID_list.py
FASTAIN=$PROJ/references/pacbio_mm10_noheader.fa
PBLIST1=$PROJ/analysis_output/pacbio_mm10_apn0.txt
PBLIST2=$PROJ/analysis_output/pacbio_mm10_apn5.txt
PBLIST3=$PROJ/analysis_output/pacbio_mm10_apn10.txt
PBOUT1=$PROJ/references/pacbio_mm10_apn0.fa
PBOUT2=$PROJ/references/pacbio_mm10_apn5.fa
PBOUT3=$PROJ/references/pacbio_mm10_apn10.fa

python2.7 $SCRIPT -i $FASTAIN -l $PBLIST1 -o $PBOUT1
python2.7 $SCRIPT -i $FASTAIN -l $PBLIST2 -o $PBOUT2
python2.7 $SCRIPT -i $FASTAIN -l $PBLIST3 -o $PBOUT3

sed -e 's/|Known|full-splice_match//g' pacbio_mm10.fa | 
sed -e 's/|Known|incomplete-splice_match//g' |
sed -e 's/|Novel|antisense//g' |
sed -e 's/|Novel|fusion//g' |
sed -e 's/|Novel|genic//g' |
sed -e 's/|Novel|genic_intron//g' |
sed -e 's/|Novel|intergenic//g' |
sed -e 's/|Novel|novel_in_catalog//g' |
sed -e 's/|Novel|novel_not_in_catalog//g' > pacbio_mm10_noheader.fa
