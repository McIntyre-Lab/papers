#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=plot
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/plot_transcriptome_MaizeWang2018_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=4
#SBATCH --array=6

module load python/3.7

## Use gffutils database (produced in the first set of EA annotations)
## 1) extract all introns
## 2) flag all exons that contain entire introns
## 3) merge flags with exon fragment information by exon ID

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
SCRIPTS=$PROJ/scripts/pacbio_transcriptome_eval
DIST=$PROJ/MaizeWang2018_transcriptome_eval/transcriptDistance/${SAMPLEID}_transcriptDistance
IND=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/${SAMPLEID}_150bp_annotations
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/check_IR
    mkdir -p ${OUTD}

## Get GFF DB
GFFDB=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.converted.gff.db

## Flag fragments involved in intron retention
python ${SCRIPTS}/flag_intron_retention.py \
    -d ${GFFDB} \
    -f ${IND}/${SAMPLEID}_exon_fragment_annotations.csv \
    -p ${SAMPLEID} \
    -o ${IND}

## Classify types of splicing in genes
##     (alt. exon, alt. donor/acceptr, intron retention)
#python $PROJ/scripts/classify_gene_splicing_02avn.py \
#    -d ${DIST}/${SAMPLEID}_pairwise_transcript_distance_reduced_sharedJunction.csv \
#    -f ${IND}/${SAMPLEID}_exon_fragment_annotations_flag_IR.csv \
#    -o ${OUTD}/${SAMPLEID}_gene_flag_splicing.csv

## Describe genes with 2 transcripts
#python $PROJ/scripts/describe_2_xcrpt_genes.py \
#    -f ${DIST}/zmtr_fsm_ism_nic_nnc_gene_flag_splicing.csv \
#    > ${DIST}/zmtr_fsm_ism_nic_nnc_describe_2_xcrpt_gene.txt
