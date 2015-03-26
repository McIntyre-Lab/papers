#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./run_identify_motifs_in_genes_up_and_down.sh 
# 
#   DESCRIPTION: This script runs identify_motifs_in_genes.pl 
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 08/22/2012 12:11:39 PM EDT
#      REVISION:  ---
#===============================================================================


PROJ=$MCLAB/arbeitman/arbeitman_fru_network
MOTIF=$PROJ/motif_analysis
GENE=$MCLAB/useful_dmel_data/flybase530/symbol2coord.csv

perl $PROJ/scripts/identify_motifs_in_genes_v6.pl -m $MOTIF/parse_mast_fru_a.csv -g $GENE -u 2000 -d 2000 >$MOTIF/fru_a_results_up_and_down.csv
perl $PROJ/scripts/identify_motifs_in_genes_v6.pl -m $MOTIF/parse_mast_fru_b.csv -g $GENE -u 2000 -d 2000 >$MOTIF/fru_b_results_up_and_down.csv
perl $PROJ/scripts/identify_motifs_in_genes_v6.pl -m $MOTIF/parse_mast_fru_c.csv -g $GENE -u 2000 -d 2000 >$MOTIF/fru_c_results_up_and_down.csv
