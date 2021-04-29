#!/bin/bash
###############################################################################
# Run Transcript Distance
# This is a bash script to describe a set of transcripts by various distance
#     measures follwing the Event Analysis annotation building step.
#
# Usage:
# Running in home directory
# bash ./run_transcriptDistance.sh ${PREFIX} ${EVENT} ${JUNC} ${FRAG} \
#                                              ${FUSION} ${OUTDIR} ${CPU}
#
# PREFIX : Prefix to append to created output files
#
# EVENT  : Path to CSV file of event_id to transcript_id to gene_id from
#          Event Analysis annotations (*_event2transcript2gene_index.csv)
#
# JUNC   : Path to CSV file of junction annotations including transcript_id
#          from Event Analysis annotations (*_annotated_junctions.csv)
#
# FRAG   : Path to CSV file of fragment annotations including transcript_id
#          from Event Analysis annotations (*_exon_fragment_annotations.csv)
#
# FUSION : Path to CSV file of fusion annotations including transcript_id
#          from Event Analysis annotations (*_fusion_annotations.csv)
#
# OUTDIR : Output directory. This will be created if it does not exist
#
# CPU    : Number of cpu to use for calculation of pairwise distances to
#          increase speed of pipeline
#
###############################################################################

## Parse command line arguments
PREFIX=${1}
EVENT=${2}
JUNC=${3}
FRAG=${4}
FUSION=${5}
OUTDIR=${6}
CPU=${7}

## Set paths
SCRIPTS=$(pwd)
if [ ! -e ${OUTDIR} ]; then mkdir -p ${OUTDIR}; fi

## Set log file - remove and replace if already exists
LOGFILE=${OUTDIR}/${PREFIX}_transcript_distance_log.txt
if [[ -f ${LOGFILE} ]]; then
	rm ${LOGFILE}
	touch ${LOGFILE}
fi

######################################
#### DESCRIBE TRANSCRIPT DISTANCE ####
######################################

## Get start time of process
start=$(date +%s.%N)
## Extract all unique pairs of individual transcript_id to gene_id
## This is the file junction variables will be merged on
echo "Extracting unique transcript_id gene_id pairs"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/extract_transcript_id_gene_id_pairs.py \
    -i ${EVENT} \
    -o ${OUTDIR}/${PREFIX}_transcript_id_2_gene_id.csv

## Get transcript-level junction variables
echo "Get transcript-level junction variables"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/get_transcript_level_junction_variables_04avn.py \
    -i ${OUTDIR}/${PREFIX}_transcript_id_2_gene_id.csv \
    -j ${JUNC} \
    -o ${OUTDIR}/${PREFIX}_transcript_level_junction.csv \
    >> ${LOGFILE}

## Get transcript-level fragment variables
echo "Get transcript-level fragment variables"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/get_transcript_level_frag_variables_03avn.py \
    -f ${FRAG} \
    -o ${OUTDIR}/${PREFIX}_transcript_level_frag.csv \
    >> ${LOGFILE}

## Get transcript-level exon region variables
echo "Get transcript-level fusion variables"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/get_transcript_level_fusion_variables_03avn.py \
    -f ${FUSION} \
    -o ${OUTDIR}/${PREFIX}_transcript_level_exonRegion.csv \
    >> ${LOGFILE}

## Merge transcript-level junction, fragment, and exon region  variables
echo "Merge transcript-level junction, fragment, and exon region variables"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/merge_transcript_level_junction_frag_exonRegion_02avn.py \
    -j ${OUTDIR}/${PREFIX}_transcript_level_junction.csv \
    -f ${OUTDIR}/${PREFIX}_transcript_level_frag.csv \
    -e ${OUTDIR}/${PREFIX}_transcript_level_exonRegion.csv \
    -o ${OUTDIR}/${PREFIX}_transcript_level_merged.csv \
    >> ${LOGFILE}


## Calculate pairwise distance measures
echo "Calculate pairwise distance measures"
/usr/bin/time -f "\tTime elapsed (mm:ss) = %E\n\tMax Mem (kb) = %M\n" \
python ${SCRIPTS}/calculate_pairwise_distance_03avn.py \
    -i ${OUTDIR}/${PREFIX}_transcript_level_merged.csv \
    -o ${OUTDIR}/${PREFIX}_pairwise_transcript_distance.csv \
    -n ${CPU} \
    >> ${LOGFILE}

## Get total time of script running
duration=$(echo "$(date +%s.%N) - $start" | bc)
hours=$(echo "$duration / 3600" | bc)
minutes=$(echo "($duration % 3600) / 60" | bc)
seconds=$(echo "($duration % 3600) % 60" | bc)
echo "
Total Runtime: $hours:$minutes:$seconds (hh:mm:ss)
"

