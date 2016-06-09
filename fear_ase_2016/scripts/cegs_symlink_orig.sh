#!/bin/bash - 
#===============================================================================
#
#         USAGE: ./cegs_symlink_orig.sh 
# 
#   DESCRIPTION: This script will crawl the original data folder structure and
#   create symlinks in a single folder so that the files can be easly accessed.
# 
#===============================================================================

PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
OUT=$PROJ/original_data/transcriptome/combined

find $PROJ/original_data/transcriptome/complete -type f -name "*.txt" -exec ln -s {} $OUT/ \;
find $PROJ/original_data/transcriptome/incomplete -type f -name "*.txt" -exec ln -s {} $OUT/ \;
find $PROJ/original_data/transcriptome/plate_k -type f -name "w67_M1*.txt" -exec ln -s {} $OUT/ \;
find $PROJ/original_data/transcriptome/plate_k -type f -name "r217_V3*.txt" -exec ln -s {} $OUT/ \;
