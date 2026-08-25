#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=extract_ujc
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/extract_ujc_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

## Description: Extract transcripts from self-map UJC GTF
module purge
ml gffread/0.12.7

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
SCRIPTS=$PROJ/scripts
REFERENCES=/blue/mcintyre/share/references

DESIGN_FILE=$PROJ/design_files/df_extract_ujc_tr_02ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

ANNO=${ARRAY[0]}
ANNO_PATH=${ARRAY[1]}
GENOME_PATH=${ARRAY[3]}

OUTD=$PROJ/annotation_id_ujc_output

# Extract transcripts

# gffread: -w=output location. -g=genome FASTA. final line = reference GTF.
gffread --w-nocds \
-w $OUTD/${ANNO}.fasta \
-g ${GENOME_PATH} \
${ANNO_PATH}

echo "Number of UJC FASTA transcripts:"
grep -c "^>" $OUTD/${ANNO}.fasta

