#!/bin/bash - 
#===============================================================================
# 
#   DESCRIPTION: This script will build the combined design files buy parsing
#   the combined folder.
# 
#===============================================================================

set -o nounset                              # Treat unset variables as an error

PROJ=/scratch/lfs/mcintyre/cegs_ase_paper
ORIG=$PROJ/original_data/transcriptome/combined

D1=$PROJ/design_files/CEGS_combined_lines_design.txt
D2=$PROJ/design_files/CEGS_combined_lines_no_tech.txt
D3=$PROJ/design_files/CEGS_combined_lines_design_r1.txt
rm $D1; rm $D2; rm $D3

for FILE in $ORIG/*
do
    NAME=$(basename $FILE .txt)
    LINE=$(echo $NAME | cut -f1 -d'_')
    MS=$(echo $NAME |sed 's/^.*_\([MV]\).*/\1/')
    REP=$(echo $NAME |sed 's/^.*_[MV]\([1-9]\+\).*/\1/')
    TECH=$(echo $NAME |sed 's/^.*_[MV][1-9]\+\.\(.*\)/\1/')
    echo "${LINE},${MS}${REP},${TECH}" >> $D1
done

cut -f1-2 -d',' $D1 | sort -u >$D2

grep -e "1$" $D1 > $D3
