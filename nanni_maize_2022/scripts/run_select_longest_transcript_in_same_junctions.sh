#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Select transcript pairs with identical junctions (prop_junction_similar == 1)
## Group by gene and set of junctions and select the longest
##     (max total nt) to represent the set

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SCRIPTS=~/mclab/SHARE/McIntyre_Lab/useful_mclab_info/scripts/transcript_distance/calculate_transcript_distance
IND=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
OUTD=$PROJ/sqanti_classification_category_subset/plot_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
    mkdir -p ${OUTD}

## Select transcripts
python $PROJ/scripts/select_longest_transcript_in_same_junctions.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance.csv \
    -d ${IND} \
    -p zmtr_fsm_ism_nic_nnc

## Plot transcript distance pairwise distances for
##     representative reference transcripts of FSM/ISM
##     and the monoexon filtered NIC/NNC after reduction of transcripts
##     with identical junctions

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction

## Split distance file by comparisons with all same junction and not
awk -F "," '{if(NR==1){ \
    for(i=1;i<=NF;i++){if($i=="prop_ER_similar"){checkCol=i; break;} \
    }print $0}else{if($checkCol==1){print $0}}}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    > ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv
awk -F "," '{if(NR==1){ \
    for(i=1;i<=NF;i++){if($i=="prop_ER_similar"){checkCol=i; break;} \
    }print $0}else{if($checkCol!=1){print $0}}}' \
    ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    > ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv


###### For the set of comparisons with all ER shared
echo "
###### For the set of comparisons with all ER shared"

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_allERshared

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_allERshared

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_allERshared


###### For the set of comparisons without all ER shared
echo "
###### For the set of comparisons without all ER shared"

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_notAllERshared

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_notAllERshared

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p zmtr_fsm_ism_nic_nnc_reduced_sharedJunction_notAllERshared

