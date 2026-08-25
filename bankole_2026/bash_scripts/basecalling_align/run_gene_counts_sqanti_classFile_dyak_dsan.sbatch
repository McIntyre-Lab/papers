#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=geneCntTO
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/geneCntTO_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --array=2-9
##SBATCH --array=2


## count genes in classification file
## need decide when to combine TRs!

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
SQANTI=${PROJ}/sqanti3QC_yak_san_dorReads_5species_dyak2_dsan1_ujc

## create output dirs
CNTS=${PROJ}/gene_transcript_counts_from_sqanti_classFile
    mkdir -p $CNTS

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
SAMPLE=${ARRAY[1]}
BC=${ARRAY[2]}
BCNUM=${ARRAY[3]}
TR=${ARRAY[4]}
TECH=${ARRAY[5]}
DATE=${ARRAY[6]}
ONT_DIR=${ARRAY[7]}

SAMPLEID=${SAMPLE}_TR${TR}

module load python/3.10

if [[ $SAMPLEID == dsan* ]]; then
    TAG=dsan1
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dsan1_ujc.gtf
else
    TAG=dyak2
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dyak2_ujc.gtf
fi


python $PROJ/scripts/gene_counts_from_classification_03amm.py \
    --classification ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt \
    --ref $REF_GTF \
    --out $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt

# add sampleID column to output
awk -v OFS=',' -v sampleid=$SAMPLEID '{if (NR==1) print "sampleID", $0; else print sampleid, $0}' $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.txt \
    > $CNTS/${SAMPLEID}_fl_rsc_${TAG}_geneCnts.csv

rm $CNTS/${SAMPLEID}_fl_rsc_dser1_geneCnts.txt
