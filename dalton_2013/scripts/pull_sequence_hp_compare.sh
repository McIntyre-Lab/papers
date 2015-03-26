#!/bin/bash - 
#===============================================================================
#
#          FILE:  pull_sequence_hp_compare.sh
# 
#         USAGE:  ./pull_sequence_hp_compare.sh 
# 
#   DESCRIPTION: For LAST pipeline test, Marty pulled a list of 81 reads that
#   were aligned by both LAST and MOSAIK. We wanted to look at these reads and
#   see if they made sense. This script pulls the sequence information and
#   count information and dumps it out for blasting.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 12/22/2011 03:48:44 PM EST
#      REVISION:  ---
#===============================================================================
WORK=/home/jfear/mclab/Fru_network/pipeline_output/data_for_checks/checking_last_accuracy
echo "read_id,alignment_quality,alignment_length,segment_of_read_that_aligned,read_length,original_read,mismatches" > $WORK/reads_in_both_mosaik_hp_3prime_and_last_hp_3prime_plus_sequence.txt;


for file in $(cat $WORK/reads_in_both_mosaik_hp_3prime_and_last_hp_3prime.txt)
do
    grep $file $WORK/read_information_uniquely_aligned_alignments_homopolymers_on_3prime_end.csv >> $WORK/reads_in_both_mosaik_hp_3prime_and_last_hp_3prime_plus_sequence.txt;
done;
