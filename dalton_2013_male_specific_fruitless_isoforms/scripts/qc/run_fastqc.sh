#!/bin/bash - 
#===============================================================================
#   DESCRIPTION:  Script to launch off fastqc qsubs
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#   INSTITUTION: University of Florida
#       CREATED: 12/02/2011 04:29:08 PM EST
#      REVISION:  ---
#===============================================================================

for DIR in /scratch/hpc/mccrory/fru_network/original_data/*
do 
    if [ -d "$DIR" ]
    then
        MYDATE=`basename $DIR`
        sed -i "s/DATE/$MYDATE/g" /scratch/hpc/jfear/Fru_network/scripts/fastqc.qsub
        qsub /scratch/hpc/jfear/Fru_network/scripts/fastqc.qsub
        sed -i "s/$MYDATE/DATE/g" /scratch/hpc/jfear/Fru_network/scripts/fastqc.qsub
    fi
done
