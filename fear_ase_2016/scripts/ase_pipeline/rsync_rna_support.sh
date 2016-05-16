#!/bin/bash
LOCALDIR=$HOME/storage/cegs_ase_paper/
HPCDIR=/scratch/lfs/mcintyre/cegs_ase_paper/ase_pipeline_output/ase_masked_aln_rna_support

rsync -av hic:$HPCDIR $LOCALDIR
