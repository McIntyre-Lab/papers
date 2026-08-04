#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=geneCntK
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/geneCntK_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --array=2-49
##SBATCH --array=2


## count genes in classification file
## need decide when to combine TRs!

## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
SQANTI=${PROJ}/sqanti3QC_ser_dorReads_2_5species_dser1_ujc

## create output dirs
CNTS=${PROJ}/gene_transcript_counts_from_sqanti_classFile
    mkdir -p $CNTS

## Get design file information (2 - 49)
DESIGN_FILE=$PROJ/design_files/design_file_sample2BC_w_origDataPath_04amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

module load python/3.10

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dser1_ujc.gtf

python $PROJ/scripts/gene_counts_from_classification_03amm.py \
    --classification ${SQANTI}/${SAMPLEID}_fl_rsc_dser1_classification.txt \
    --ref $REF_GTF \
    --out $CNTS/${SAMPLEID}_fl_rsc_dser1_geneCnts.txt

# add sampleID column to output
awk -v OFS=',' -v sampleid=$SAMPLEID '{if (NR==1) print "sampleID", $0; else print sampleid, $0}' $CNTS/${SAMPLEID}_fl_rsc_dser1_geneCnts.txt \
    > $CNTS/${SAMPLEID}_fl_rsc_dser1_geneCnts.csv

rm $CNTS/${SAMPLEID}_fl_rsc_dser1_geneCnts.txt


