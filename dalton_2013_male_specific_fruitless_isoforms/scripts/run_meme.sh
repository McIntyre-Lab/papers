#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./run_meme.sh 
# 
#   DESCRIPTION: This script runs Meme Suite's MEME program on sequences
#   simulated from Position Weight Matrices.
# 
#  REQUIREMENTS: meme.bin
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 08/20/2012 02:47:55 PM EDT
#      REVISION:  ---
#===============================================================================

MEME=/home/jfear/opt/meme/bin/meme.bin
PROJ=$MCLAB/arbeitman_fru_network
ORIG=$PROJ/data

$MEME $ORIG/fru_a_simulated_sequences_jmf.fa -o $PROJ/motif_analysis/meme_fru_a
$MEME $ORIG/fru_b_simulated_sequences_jmf.fa -o $PROJ/motif_analysis/meme_fru_b
$MEME $ORIG/fru_c_simulated_sequences_jmf.fa -o $PROJ/motif_analysis/meme_fru_c
