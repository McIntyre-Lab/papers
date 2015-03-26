#!/bin/bash - 
#===============================================================================
#
#          FILE:  make_barcode_table.sh
# 
#         USAGE:  ./make_barcode_table.sh 
# 
#   DESCRIPTION:  Strip barcode information out of the lanes for later use.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 11/12/2011 09:23:05 AM EST
#      REVISION:  ---
#===============================================================================

WORK='/scratch/hpc/jfear/fru_network/Data';
echo "lane,seq_num,barcode" > /scratch/hpc/jfear/fru_network/barcode_table.csv;

for i in {1..8}
do
    for j in {1..4}
    do
        if [ -e $WORK/Lane$i/s_${i}_00${j}_sequence.txt ]
        then
            BAR=`head -n 1 $WORK/Lane$i/s_${i}_00${j}_sequence.txt | sed 's/.*\#\(.*\)\/1/\1/g'` ;      
            echo "$i,00$j,$BAR" >> /scratch/hpc/jfear/fru_network/barcode_table.csv;
        fi;
    done;
done;
