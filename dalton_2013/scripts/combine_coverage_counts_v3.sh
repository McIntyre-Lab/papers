#!/bin/bash - 
#===============================================================================
#
#   DESCRIPTION:  This script will combine all of the coverage counts into one
#   large file. 
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 11/12/2011 04:02:35 PM EST
#      REVISION: 12/21/2011
#                 - Edited the script for use with a coverage counts from a new
#                   alignment pipeline "coverage_on_fusions"
#                 - Also in the initial coverage counts files chromosome and
#                   position information was included. The previous version
#                   removed this information with a sed statment. These coverage
#                   counts did not have this so this statement was removed.
#                 - Added header information to the output file
#                12/23/2011
#                 - Changed line 27 so that only the files I want are included
#===============================================================================

WORK="$MCLAB/Fru_network/pipeline_output/coverage_on_fusions";

cd $WORK;

echo "fusion_id,sample_id,mapped_reads,reads_in_exon,coverage_in_exon,exon_length,apn,rpkm" > $WORK/all_coverage_counts.csv;
for file in coverage_on_fusions_2011-05-03_1_ACTTGA.csv coverage_on_fusions_2011-05-03_1_ATCACG.csv coverage_on_fusions_2011-05-03_1_CAGATC.csv coverage_on_fusions_2011-05-03_1_GATCAG.csv coverage_on_fusions_2011-05-03_1_TAGCTT.csv coverage_on_fusions_2011-05-03_2_ACTTGA.csv coverage_on_fusions_2011-05-03_2_ATCACG.csv coverage_on_fusions_2011-05-03_2_GATCAG.csv coverage_on_fusions_2011-05-03_2_GCCAAT.csv coverage_on_fusions_2011-05-03_2_TAGCTT.csv coverage_on_fusions_2011-05-03_3_ACAGTG.csv coverage_on_fusions_2011-05-03_3_ACTTGA.csv coverage_on_fusions_2011-05-03_3_GATCAG.csv coverage_on_fusions_2011-05-03_3_TAGCTT.csv coverage_on_fusions_2011-05-03_3_TGACCA.csv coverage_on_fusions_2011-05-03_5_ACTTGA.csv coverage_on_fusions_2011-05-03_5_CGATGT.csv coverage_on_fusions_2011-05-03_5_GATCAG.csv coverage_on_fusions_2011-05-03_5_TGACCA.csv coverage_on_fusions_2011-05-03_5_TTAGGC.csv coverage_on_fusions_2011-05-03_6_ACAGTG.csv coverage_on_fusions_2011-05-03_6_ACTTGA.csv coverage_on_fusions_2011-05-03_6_CGATGT.csv coverage_on_fusions_2011-05-03_6_TAGCTT.csv coverage_on_fusions_2011-05-03_6_TTAGGC.csv coverage_on_fusions_2011-05-03_7_TAGCTT.csv coverage_on_fusions_2011-07-05_1_ATCACG.csv coverage_on_fusions_2011-07-05_1_CAGATC.csv coverage_on_fusions_2011-07-05_1_GCCAAT.csv coverage_on_fusions_2011-07-05_1_TGACCA.csv coverage_on_fusions_2011-07-05_1_TTAGGC.csv coverage_on_fusions_2011-07-05_2_ACAGTG.csv coverage_on_fusions_2011-07-05_2_ACTTGA.csv coverage_on_fusions_2011-07-05_2_ATCACG.csv coverage_on_fusions_2011-07-05_2_CGATGT.csv coverage_on_fusions_2011-07-05_2_TGACCA.csv coverage_on_fusions_2011-07-05_3_ACAGTG.csv coverage_on_fusions_2011-07-05_3_CGATGT.csv coverage_on_fusions_2011-07-05_3_GATCAG.csv coverage_on_fusions_2011-07-05_3_GCCAAT.csv coverage_on_fusions_2011-07-05_3_TTAGGC.csv coverage_on_fusions_2011-07-05_5_ATCACG.csv coverage_on_fusions_2011-07-05_5_CAGATC.csv coverage_on_fusions_2011-07-05_5_CGATGT.csv coverage_on_fusions_2011-07-05_5_TAGCTT.csv coverage_on_fusions_2011-07-05_5_TTAGGC.csv coverage_on_fusions_2011-07-05_6_ACTTGA.csv coverage_on_fusions_2011-07-05_6_ATCACG.csv coverage_on_fusions_2011-07-05_6_CGATGT.csv coverage_on_fusions_2011-07-05_6_TAGCTT.csv coverage_on_fusions_2011-07-05_6_TGACCA.csv coverage_on_fusions_2011-07-05_7_ACAGTG.csv coverage_on_fusions_2011-07-05_7_ACTTGA.csv coverage_on_fusions_2011-07-05_7_CAGATC.csv coverage_on_fusions_2011-07-05_7_GATCAG.csv coverage_on_fusions_2011-07-05_8_ACAGTG.csv coverage_on_fusions_2011-07-05_8_ATCACG.csv coverage_on_fusions_2011-07-05_8_CGATGT.csv coverage_on_fusions_2011-07-05_8_GATCAG.csv coverage_on_fusions_2011-07-05_8_GCCAAT.csv 
do 
    if [[ $file != *unknown* ]] 
    then 
        cat $file >> $WORK/all_coverage_counts.csv 
        echo "$file done";
    fi 
done;


