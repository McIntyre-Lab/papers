#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./run_parse_mast_xml.sh 
# 
#   DESCRIPTION: This script runs the perl scrip that parses the mast xml
#   output.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 08/21/2012 03:52:06 PM EDT
#      REVISION:  ---
#===============================================================================

PROJ=$MCLAB/arbeitman_fru_network
SCRIPT=$PROJ/scripts/parse_mast_xml_output_v2.pl
OUTPUT=$PROJ/motif_analysis

$SCRIPT $OUTPUT/mast_fru_a/mast.xml > $OUTPUT/parse_mast_fru_a.csv
$SCRIPT $OUTPUT/mast_fru_b/mast.xml > $OUTPUT/parse_mast_fru_b.csv
$SCRIPT $OUTPUT/mast_fru_c/mast.xml > $OUTPUT/parse_mast_fru_c.csv
