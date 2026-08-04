#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=geneCntR
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/rmg_lmm_dros_data/scripts/SLURM_LOGS/geneCntR_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
##SBATCH --array=2-109
#SBATCH --array=2-6


## count genes in classification file
## need decide when to combine TRs!

## input dirs
PROJ=/blue/mcintyre/share/rmg_lmm_dros_data
SQANTI=${PROJ}/sqanti3QC_mel_sim_dorReads_2_5species_dmel6_dsimW_ujc/

## create output dirs
CNTS=${PROJ}/gene_transcript_counts_from_sqanti_classFile_test
    mkdir -p $CNTS

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_bc_design_w_origDataPath_04amm.csv
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

if [[ $SAMPLEID == mel* ]]; then
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf
    TAG=dmel6
else
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dsimW_ujc.gtf
    TAG=dsimW
fi

python $PROJ/scripts/gene_counts_from_classification_02ndk.py \
    --classification ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt \
    --ref $REF_GTF \
    --out $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt


# add sampleID column to output
awk -v OFS=',' -v sampleid=$SAMPLEID '{if (NR==1) print "sampleID", $0; else print sampleid, $0}' $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt \
    > $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.csv

rm $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt


