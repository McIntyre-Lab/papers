#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Use gffutils database (produced in the first set of EA annotations)
## 1) extract all introns
## 2) flag all exons that contain entire introns
## 3) merge flags with exon fragment information by exon ID

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=$PROJ/scripts
IND=$PROJ/sqanti_classification_category_subset/EA_annotations

## Get GFF DB
GFFDB=$PROJ/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.EA_converted.gff.db

## Flag fragments involved in intron retention
python ${SCRIPTS}/flag_intron_retention.py \
    -d ${GFFDB} \
    -f ${IND}/zmtr_fsm_ism_nic_nnc_exon_fragment_annotations.csv \
    -p zmtr_fsm_ism_nic_nnc \
    -o ${IND}
