#!/bin/bash - 
#===============================================================================
#
#   DESCRIPTION:  I did not realize that when Marty added sample ID to these
#   file he kept the number format for the 08-25-2010 series. So I need to edit
#   this so that instead of 001-004 they have the barcode IDs. I can simply
#   take this from the file name now.
# 
#  REQUIREMENTS:  barcod_table.csv
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 11/12/2011 10:22:07 AM EST
#      REVISION:  ---
#===============================================================================

WORK='/share/mclab/Fru_network/pipeline_output';

# you can use either
# $WORK/coverage
# $WORK/coverage_spikes

cd $WORK/coverage_spikes;
for file in coverage_2010-08-25_*
do 
    FILENAME=`basename $file .csv`;
    LANE=`echo $FILENAME | sed 's/coverage_2010-08-25_\([1-8]\)_.*/\1/g' `
    BAR=`echo $FILENAME | sed 's/coverage_2010-08-25_[1-8]_\(.*\)/\1/g' `
    sed -i "s/\(^.*,2010-08-25_${LANE}_\)00[1-9]\(,.*$\)/\1${BAR}\2/g" $file;
    echo "$FILENAME done";
done;
