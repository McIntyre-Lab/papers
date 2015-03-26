#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./ 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 10/03/2012 10:46:14 AM EDT
#      REVISION:  ---
#===============================================================================

PROJ=$MCLAB/arbeitman_fru_network
SCRIPTS=$PROJ/scripts
MEANS=$PROJ/data/for_wiggles

FILE=$MEANS/AH_CS.csv
NAME=`basename $FILE .csv`
SYM2CORD=$MCLAB/useful_dmel_data/flybase530/symbol2coord.csv

GENE=Ace

$SCRIPTS/parse_gene_pileup.pl -c $SYM2CORD -f $FILE -g $GENE >$HOME/tmp/$NAME.$GENE.csv

