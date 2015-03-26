#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./create_simulated_motifs.sh 
# 
#   DESCRIPTION: This script creates the simulated motifs from a Position
#   Weight Matrix (PWM)
# 
#  REQUIREMENTS: generate_simulated_meme_sequences_v2_jmf.pl
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 08/20/2012 02:42:13 PM EDT
#      REVISION:  ---
#===============================================================================

# PWMs to use for simulation
FRUA=$MCLAB/arbeitman_fru_network/original_data/fru_a_weight_matrix.tsv
FRUB=$MCLAB/arbeitman_fru_network/original_data/fru_b_weight_matrix.tsv
FRUC=$MCLAB/arbeitman_fru_network/original_data/fru_c_weight_matrix.tsv

perl $MCLAB/arbeitman_fru_network/scripts/generate_simulated_meme_sequences_v3_jmf.pl $FRUA > $MCLAB/arbeitman_fru_network/data/fru_a_simulated_sequences_jmf.fa
perl $MCLAB/arbeitman_fru_network/scripts/generate_simulated_meme_sequences_v3_jmf.pl $FRUB > $MCLAB/arbeitman_fru_network/data/fru_b_simulated_sequences_jmf.fa
perl $MCLAB/arbeitman_fru_network/scripts/generate_simulated_meme_sequences_v3_jmf.pl $FRUC > $MCLAB/arbeitman_fru_network/data/fru_c_simulated_sequences_jmf.fa
