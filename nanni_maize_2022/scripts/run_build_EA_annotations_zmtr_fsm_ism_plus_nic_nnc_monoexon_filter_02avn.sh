#!/bin/sh

## Example shell script to build Event Analysis annotations

# Variables to set

## Path to Event Analysis install
EVENTDIR=~/event_analysis/events/event_analysis

## Label to prefix to annotation files (e.g. mm10_refseq, dmel617, hg38_ens)
PREFIX=zmtr_fsm_ism_nic_nnc

## Path to formatted GFF3 file (should in a FlyBase GFF3 format)
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
GFF=$PROJ/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter.EA_converted.gff

## Path to genome FASTA
FASTA=~/mclab/SHARE/McIntyre_Lab/useful_maize_info/Zea_mays.B73_RefGen_v4.dna.toplevel.fa


## Size (in nt) of reads or each read pair (set to maximum if you have reads of 
## various lengths, i.e. from adapter trimming)
READSIZE=150

## Output directory. If it does not exist, the annotation build script with create it
OUTDIR=/TB14/TB14/roz_maize_ozone_EA

# Build Event Analysis annotations
cd ${EVENTDIR}
sh ./run_buildAnnotations.sh ${PREFIX} ${GFF} ${FASTA} ${READSIZE} ${OUTDIR}

rsync -urlvz ${OUTDIR}/ $PROJ/sqanti_classification_category_subset/EA_annotations/.

rm -r ${OUTDIR}
