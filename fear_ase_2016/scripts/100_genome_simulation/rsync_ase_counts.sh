#!/bin/bash

LOCALDIR=/home/jfear/mclab/cegs_ase_paper/pipeline_output/ase_counts_fb557_100_genome_simulation
REMOTEDIR=/scratch/lfs/mcintyre/cegs_ase_paper/100_genome_simulation/ase_counts_fb557

if [ ! -e $LOCALDIR ]; then mkdir -p $LOCALDIR; fi;

rsync -rultv hic:$REMOTEDIR/* $LOCALDIR/
