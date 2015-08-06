#!/bin/bash - 
#===============================================================================
#
#   DESCRIPTION: Sync the simulation SAS DATA to the local machine for
#   summarizing the results.
# 
#        AUTHOR: Justin Fear (JMF), jfear@ufl.edu
#  ORGANIZATION: University of Florida
#       CREATED: 01/15/2015 08:11:02 AM EST
#===============================================================================

set -o nounset                              # Treat unset variables as an error

LOCALDIR=$HOME/tmp
REMOTEDIR=/scratch/lfs/mcintyre/cegs_sem_sd_paper/cegs_adding_genes_simulation

NUMGENES=800
NUMSIMS=10


for I in `seq 1 $NUMSIMS`
do
    if [ ! -e $LOCALDIR/${NUMGENES}_gene_sim/$I ]; then mkdir -p $LOCALDIR/${NUMGENES}_gene_sim/$I; fi

    rsync -av hic:$REMOTEDIR/${NUMGENES}_gene_sim/$I/sas_data $LOCALDIR/${NUMGENES}_gene_sim/$I/

done
