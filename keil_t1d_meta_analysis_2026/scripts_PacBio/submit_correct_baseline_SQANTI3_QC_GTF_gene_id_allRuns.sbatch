#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=correct_gene_id
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/correct_gene_id_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --array=2-13

## Correct the gene_id values of the GTF output from SQANTI3 QC
##     to match the gene_id values associated with the
##     transcripts in the classification file

module load python/3.6

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
SCRIPTS=$PROJ/scripts
IND=$PROJ/sqanti3_qc_allRuns

## design has header!!
DESIGN_FILE=$PROJ/design_files/design_sampleID_2_barcode.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
BC=${ARRAY[2]}

date

echo "
*** Correct gene_id values ***

  inputs:  ${IND}/${SAMPLEID}_${RUN}_classification.tx
           {IND}/${SAMPLEID}_${RUN}_corrected.gtf
"
for RUN in first second third fourth fifth
do
   
    ## Get associated gene values for each isoform from classification files
    ## Correct the gene_id values of the GTF output from SQANTI3 QC to
    ##     match the gene_id values associated with the transcripts in
    ##     the classification file

    python ${SCRIPTS}/correct_SQANTI_QC_GTF_gene_id.py \
        -c ${IND}/${SAMPLEID}_${RUN}_classification.txt \
        -g ${IND}/${SAMPLEID}_${RUN}_corrected.gtf \
        -o ${IND}/${SAMPLEID}_${RUN}_corrected_associated_gene.gtf

done
