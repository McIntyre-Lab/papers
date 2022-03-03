#!/bin/bash

## Convert filtered and corrected GTF to Flybase GFF3 for creating annotations
## Run under conda environment with EA v1.0.16 installed
## JRBN on ufgi-120c-w05 run: `conda activate ea1016` FIRST!


## getting database locked errors, so going to try syncing input GTF to my local do what I need to there
## then move final output back to the share

PROJ=$MCLAB/maize_ozone_FINAL

GTFIN=$PROJ/2018/PacBio/sqanti_post_filter_b73/sqanti_b73_filtered_corrected.gtf

GFFOUT=$PROJ/2018/PacBio/sqanti_post_filter_b73/sqanti_b73_filtered_corrected.EA_converted.gff

EVENTDIR=$MCLAB/jrbnewman/event_analysis

SCRIPTS=$EVENTDIR/event_analysis/src

ROZ=$HOME/tmp
mkdir -p $ROZ

rsync -rtDcv $GTFIN $ROZ

GTFIN2=$ROZ/sqanti_b73_filtered_corrected.gtf
GFFOUT2=$ROZ/sqanti_b73_filtered_corrected.EA_converted.gff

# Checking if a gff.db file exists, and if not then create one
echo "Checking if user-supplied GTF/GFF file has a pre-generated database file"
    if [ ! -e ${GTFIN2}.db ]
    then
    echo "No database file found! Generating..."
    python $SCRIPTS/make_gff_db.py --gff $GTFIN2
    echo "Database generated!"
    else
    echo "Database file found! Skipping generation."
    fi

# Convert user-supplied GTF/GFF file into GFF3 format
echo "Converting user-supplied GTF/GFF to GFF3 format"

python $SCRIPTS/convertGTF2GFF3.py --input $GTFIN2 --output ${GFFOUT2}

sort -k1,1 -k4n -k5n ${GFFOUT2} > ${ROZ}/temp.gff
mv ${ROZ}/temp.gff ${GFFOUT2}

python $SCRIPTS/make_gff_db.py --gff ${GFFOUT2}

rsync -rtDcv ${GFFOUT2}* $PROJ/2018/PacBio/sqanti_post_filter_b73/.
