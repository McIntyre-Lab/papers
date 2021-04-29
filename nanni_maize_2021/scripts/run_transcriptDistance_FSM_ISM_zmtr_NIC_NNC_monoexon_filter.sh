#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Shell script to describe a set of transcripts by various distance measures

#####################################################################
##################### Variables for user to set #####################
#####################################################################

## Path to transcript distance scripts
SCRIPTS=~/mclab/SHARE/McIntyre_Lab/useful_mclab_info/scripts/transcript_distance/calculate_transcript_distance

## Label to prefix to otuput files (e.g. mm10_refseq, dmel617, hg38_ens)
PREFIX=zmtr_fsm_ism_nic_nnc

## Output directory. If it does not exist, the script will create it
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
OUTDIR=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter

## Path to Event Analysis annotation output
EVENTPATH=$PROJ/sqanti_classification_category_subset/EA_annotations

## Prefix used in Event Analysis annotation files (e.g. mm10_refseq, dmel617, hg38_ens)
EVENTPREFIX=zmtr_fsm_ism_nic_nnc

## Number of CPU to use for paralellization of distance calculations
CPU=8

#####################################################################
############# Code below does not need changed by user ##############
#####################################################################

## Path to CSV file of event_id to transcript_id to gene_id from
##     Event Analysis annotations (*_event2transcript2gene_index.csv)
EVENT=${EVENTPATH}/${EVENTPREFIX}_event2transcript2gene_index.csv

## Path to CSV file of junction annotations including transcript_id
##     from Event Analysis annotations (*_annotated_junctions.csv)
JUNC=${EVENTPATH}/${EVENTPREFIX}_annotated_junctions.csv

## Path to CSV file of fragment annotations including transcript_id
##     from Event Analysis annotations (*_exon_fragment_annotations.csv)
FRAG=${EVENTPATH}/${EVENTPREFIX}_exon_fragment_annotations.csv

## Path to CSV file of fusion annotations including transcript_id
##     from Event Analysis annotations (*_fusion_annotations.csv)
FUSION=${EVENTPATH}/${EVENTPREFIX}_fusion_annotations.csv

## Collapse transcriptome by junctions
cd ${SCRIPTS}
bash ./run_transcriptDistance.sh ${PREFIX} ${EVENT} ${JUNC} ${FRAG} ${FUSION} ${OUTDIR} ${CPU}

