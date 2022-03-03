#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=select
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/select_longest_MaizeWang2018_%A_%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=8
#SBATCH --array=1-6

module load python/3.7

## Select transcript pairs with identical junctions (prop_junction_similar == 1)
## Group by gene and set of junctions and select the longest
##     (max total nt) to represent the set

## Get info from design file
PROJ=/blue/mcintyre/share/maize_ainsworth
DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Set directories
IND=$PROJ/MaizeWang2018_transcriptome_eval/transcriptDistance/${SAMPLEID}_transcriptDistance
SCRIPTS=$PROJ/scripts/transcript_distance/calculate_transcript_distance
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome
    mkdir -p ${OUTD}

## Select transcripts
python $PROJ/scripts/pacbio_transcriptome_eval/select_longest_transcript_in_same_junctions.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance.csv \
    -d ${IND} \
    -p ${SAMPLEID} \
    > ${OUTD}/${SAMPLEID}_consol_counts.txt

## Plot transcript distance pairwise distances for
##     representative reference transcripts of FSM/ISM
##     and the monoexon filtered NIC/NNC after reduction of transcripts
##     with identical junctions

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction

## Split distance file by comparisons with all same junction and not
awk -F "," '{if(NR==1){ \
    for(i=1;i<=NF;i++){if($i=="prop_ER_similar"){checkCol=i; break;} \
    }print $0}else{if($checkCol==1){print $0}}}' \
    ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
    > ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv
awk -F "," '{if(NR==1){ \
    for(i=1;i<=NF;i++){if($i=="prop_ER_similar"){checkCol=i; break;} \
    }print $0}else{if($checkCol!=1){print $0}}}' \
    ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
    > ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv


###### For the set of comparisons with all ER shared
echo "
###### For the set of comparisons with all ER shared"

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_allERshared

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_allERshared

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_allERshared.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_allERshared


###### For the set of comparisons without all ER shared
echo "
###### For the set of comparisons without all ER shared"

## Plot histograms of distance similarities for genes with 2 transcripts
echo "Plot histograms of distance similarities for genes with 2 transcripts"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 2 \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_notAllERshared

## Plot histograms of minimum distance similarities for all genes
echo "Plot histograms of minimum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 0 \
    -m min \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_notAllERshared

## Plot histograms of maximum distance similarities for all genes
echo "Plot histograms of maximum distance similarities for all genes"
python ${SCRIPTS}/plot_transcriptDistance.py \
    -i ${IND}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction_notAllERshared.csv \
    -n 0 \
    -m max \
    -c similarity \
    -d ${OUTD} \
    -p ${SAMPLEID}_reduced_sharedJunction_notAllERshared
