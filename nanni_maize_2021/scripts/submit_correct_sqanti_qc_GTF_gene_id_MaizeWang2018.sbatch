#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=correct_gene_id
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/correct_gene_id_MaizeWang2018_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

## Correct the gene_id values of the GTF output from SQANTI3 QC
##     to match the gene_id values associated with the
##     transcripts in the classification file

module load python/3.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
SCRIPTS=$PROJ/scripts/pacbio_transcriptome_eval
IND=$PROJ/MaizeWang2018_transcriptome_eval/sqanti_qc

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

date

echo "
*** Correct gene_id values ***

Inputs: ${IND}/${SAMPLEID}_classification.txt
        ${IND}/${SAMPLEID}.collapsed_corrected.gtf
"
## Get associated gene values for each isoform from classification files
## Correct the gene_id values of the GTF output from SQANTI3 QC to
##     match the gene_id values associated with the transcripts in
##     the classification file

python ${SCRIPTS}/correct_SQANTI_QC_GTF_gene_id_02avn.py \
    -c ${IND}/${SAMPLEID}_classification.txt \
    -g ${IND}/${SAMPLEID}.collapsed_corrected.gtf \
    -o ${IND}/${SAMPLEID}_corrected_associated_gene.gtf

date
