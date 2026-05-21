#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=fiveSpecies_02
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/fiveSpecies_02_%A-%a.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account=mcintyre
#SBATCH --array=2-6

module purge

## Description: Run script to create 5 species annotation, do GFFCompare to add geneIDs.

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
DESIGN_FILE=$PROJ/design_files/df_create_fiveSpecies_01ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

GENOME=${ARRAY[0]}
REF_GTF=${ARRAY[1]}
INFILES=${ARRAY[2]}

OUTD=$PROJ/fiveSpecies_annotations
	mkdir -p $OUTD

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities
NAME=fiveSpecies_2_${GENOME}_ujc

# Combine multiple UJC GTFs on the same genome, accounts for differences between GTFs at the very beginning and end of the same transcript
python $UTIL/combine_multi_ujc_gtf_on_one_genome_01ksb.py \
	-ig $INFILES \
	-o ${OUTD}/${NAME}_roz.gtf

module load gffcompare/0.12.2

# GFFCOMPARE TO GET REFERENCE GENE IDs IN OUTPUT GTF!!

gffcompare -r $REF_GTF -o $OUTD/${NAME}_roz $OUTD/${NAME}_roz.gtf

python $UTIL/correct_gffcompare_GTF_gene_id.py \
-a $OUTD/${NAME}_roz.annotated.gtf \
-g $OUTD/${NAME}_roz.gtf \
-o $OUTD/${NAME}.gtf \
--key ${OUTD}/${NAME}_gene_key.csv \
--change-XLOC GFF5

rm $OUTD/*$GENOME*roz*
