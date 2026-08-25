#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=extract_anno
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/extract_anno_02_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

## Description: Extract transcripts from reference annotation

module purge
ml gffread/0.12.7
ml samtools/1.18

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
SCRIPTS=$PROJ/scripts
REFERENCES=/blue/mcintyre/share/references

DESIGN_FILE=$PROJ/design_files/df_extract_anno_tr_03ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

ANNO=${ARRAY[0]}
ANNO_PATH=${ARRAY[1]}
GENOME=${ARRAY[2]}
GENOME_PATH=${ARRAY[3]}
REF_TR=${ARRAY[4]}

# Extract transcripts using GFFRead

# gffread: -w=output location. -g=genome FASTA. final line = reference GTF.
gffread --w-nocds \
-w ${REF_TR}.fasta \
-g ${GENOME_PATH} \
${ANNO_PATH}

# Create fasta index
samtools faidx ${REF_TR}.fasta

# Print number of transcript
echo "Number of FASTA transcripts:"
grep -c "^>" ${REF_TR}.fasta


