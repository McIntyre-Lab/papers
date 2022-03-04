#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=htseq_count
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/htseq_count_genes_loci_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --array=1-11

## Get read counts over gene regions of trimmed CCS using HTSeq

module purge
module load htseq/0.11.2

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/check_isoseq3_lost_genes/minimap2
OUTD=$PROJ/check_isoseq3_lost_genes/htseq_gene_counts
    mkdir -p ${OUTD}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

## GTF for quantification of 12604 assembled transcriptome
GTF=$PROJ/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf

## Quantify genes with htseq
htseq-count \
    ${IND}/${SAMPLE_NAME}_minimap2_b73.sam \
    ${GTF} \
    -f sam \
    -s no \
    -m union \
    --secondary-alignment=ignore \
    --supplementary-alignments=ignore \
    > ${OUTD}/${SAMPLE_NAME}_ccs_gene_htseq_count.tsv
