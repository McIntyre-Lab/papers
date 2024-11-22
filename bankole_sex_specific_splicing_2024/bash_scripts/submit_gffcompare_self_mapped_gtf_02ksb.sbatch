#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gffc_map2one
#SBATCH --output=/blue/mcintyre/share/sex_specific_splicing/scripts/SLURM_LOGS/gffc_map2one_%A-%a.out
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --account=mcintyre
#SBATCH --array=2-7

## Description: GFFCompare on self-mapped reference. Output any genes involved in gffcompare mismatches (transcript is assigned a different gene than it originally had pre-mapping).

module purge

ml gffcompare/0.12.2
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

### Set Directories
PROJ=/blue/mcintyre/share/sex_specific_splicing
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

DESIGN_FILE=$PROJ/design_files/df_gffcomp_and_id_ujc_03ksb.csv

DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

OUTNAME=${ARRAY[0]}
IND=${ARRAY[1]}
ANNO_PATH=${ARRAY[2]}
PC_REF=${ARRAY[3]}

echo "Running $OUTNAME:"

# GFFCOMPARE TO GET REFERENCE GENE IDS IN OUTPUT GTF
gffcompare -r $ANNO_PATH -o $IND/${OUTNAME}_roz $IND/${OUTNAME}_noGeneID.gtf

# AVN Utility to use GFFCOMPARE output to create the GTF with geneIDs and create a gene key
python $UTIL/correct_gffcompare_GTF_gene_id.py \
-a $IND/${OUTNAME}_roz.annotated.gtf \
-g $IND/${OUTNAME}_noGeneID.gtf \
-o $IND/${OUTNAME}.gtf \
--key ${IND}/${OUTNAME}_gffcompare_gene_key.csv \
--change-XLOC GFF

# Get list of ambig genes and list of transcript with flag for misassigned genes

python $UTIL/compare_gene_assignments_post_gffcompare_01ksb.py \
-a  $ANNO_PATH \
-pc $PC_REF \
-k  ${IND}/${OUTNAME}_gffcompare_gene_key.csv \
-p  $OUTNAME \
-o  $IND

rm $IND/$OUTNAME*roz*
rm $IND/${OUTNAME}_noGeneID.gtf
