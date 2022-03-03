#!/bin/bash

## Convert filtered and corrected GTF to Flybase GFF3 for creating annotations
## First run by JRBN under conda environment with EA v1.0.16 installed
##     JRBN on ufgi-120c-w05 run: `conda activate ea1016` FIRST!
## Rerun by AVN (2/21/21) under conda environment event_analysis_updated

## getting database locked errors, so going to try syncing input GTF to my local do what I need to there
## then move final output back to the share

PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio

GTFIN=$PROJ/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.gtf

GFFOUT=$PROJ/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.EA_converted.gff

EVENTDIR=~/event_analysis/events/event_analysis

SCRIPTS=$EVENTDIR/src

ROZ=/TB14/TB14/roz_maize_ozone_EA
mkdir -p $ROZ

rsync -rtDcv $GTFIN $ROZ


GTFIN2=$ROZ/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.gtf
GFFOUT2=$ROZ/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.EA_converted.gff

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

rsync -rtDcv ${GFFOUT2}* $PROJ/sqanti_classification_category_subset/.

rm -r ${ROZ}
