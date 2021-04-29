#!/bin/sh

## Set paths
PROJ=/blue/mcintyre/share/maize_ainsworth
DATA=/blue/mcintyre/share/transcript_distance/MaizeWang2018
DESIGN=$PROJ/design_files
INFILE=${DESIGN}/MaizeWang2018_subread_best_guess_bas_h5.csv
OUTFILE=${DESIGN}/df_MaizeWang2018_sorghum_leaf_subread.csv
NOHEADER=${DESIGN}/df_MaizeWang2018_sorghum_leaf_subread_noHeader.csv

## Extract from best guess file the pool 1 sburead fullpath and filename
awk -F "," '$4==2||NR==1{print $1","$2}' ${INFILE} > ${OUTFILE}
awk -F "," '$4==2{print $1","$2}' ${INFILE} > ${NOHEADER}

## Count
echo "
$(wc -l ${NOHEADER} | awk '{print $1}') *.bas.h5 files in design file
$(for MOVIE in $(cut -d "," -f 2 ${NOHEADER}); do find ${DATA} -name *${MOVIE}*.subreads.fastq; done | wc -l) corresponding *.subreads.fastq files
"
