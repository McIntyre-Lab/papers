#!/bin/bash - 
#===============================================================================
#
#   DESCRIPTION:  This script will combine all of the coverage counts spikes
#   into one large file. Also in the coverage count file there is chromosome
#   information that I do not want right now, so I will split this out to a new
#   file.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 11/12/2011 04:02:35 PM EST
#      REVISION:  ---
#===============================================================================

WORK='/share/mclab/Fru_network/pipeline_output/coverage_spikes';

cd $WORK;

echo "ercc_id,sample_id,mapped_reads,reads_in_exon,coverage_in_exon,exon_length,apn,rpkm" > $WORK/all_coverage_spikes.csv;
for file in *
do 
    if [[ $file != *unknown* ]] 
    then 
        cat $file >> $WORK/all_coverage_spikes.csv 
        echo "$file done";
    fi 
done;



