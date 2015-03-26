#!/bin/bash - 
#===============================================================================
#
#          FILE:  convert_number2barcode.sh
# 
#   DESCRIPTION:  For the 08/25/2010 run instead of using barcodes in the file
#   name, number 001-004 where used. Just so everything is named the same
#   barcode information was stripped out of the original data and now the
#   coverage count file are going to be renamed to have this information
# 
#  REQUIREMENTS:  barcod_table.csv
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 11/12/2011 10:22:07 AM EST
#      REVISION:  ---
#===============================================================================

WORK='/share/mclab/Fru_network/pipeline_output';

# you can use either:
# $WORK/coverage
# $WORK/coverage_spikes

cd $WORK/coverage_spikes;
for file in coverage_2010-08-25_*
do 
    FILENAME=`basename $file .csv`;
    LANE=`echo $FILENAME | sed 's/coverage_2010-08-25_\([1-8]\)_00[1-4]/\1/g' `
    NUM=`echo $FILENAME | sed 's/coverage_2010-08-25_[1-8]_\(00[1-4]\)/\1/g' `
    BAR=`grep "$LANE,$NUM" $WORK/barcode_table.csv | cut -d "," -f 3`
    mv $file coverage_2010-08-25_${LANE}_${BAR}.csv
done;
