#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=prep
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/rnaseq_prep_4_tappas_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Using on/off flags (TPM>5 in at least 50% of reps)
## Select transcripts that are detected in at least one sample
## Select columns that apply to each genotype
## Exclude "bad" libraries: B73_P1_C7_Ele B73_P4_C1_Amb B73_P4_C6_Amb

module purge
module load python/3.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INDIR=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc


## Reformat expression matrix and make design file for each genotype
## For expected_count matrix
#python $PROJ/scripts/pacbio_transcriptome_eval/prep_rnaseq_4_tappas.py \
#    -t ${INDIR}/flag_on_off_rsem_expression_isoforms_TPM_5_full_table.tsv \
#    -c ${INDIR}/combined_isoforms_expression_matrix.expected_count.txt \
#    -e B73_P1_C7_Ele -e B73_P4_C1_Amb -e B73_P4_C6_Amb \
#    -v expected_count \
#    -o ${INDIR}
## For TPM matrix
python $PROJ/scripts/pacbio_transcriptome_eval/prep_rnaseq_4_tappas_02avn.py \
    -t ${INDIR}/flag_on_off_rsem_expression_isoforms_TPM_5_full_table.tsv \
    -c ${INDIR}/combined_isoforms_expression_matrix.expected_count.txt \
    -e B73_P1_C7_Ele -e B73_P4_C1_Amb -e B73_P4_C6_Amb \
    -v TPM \
    -o ${INDIR}
