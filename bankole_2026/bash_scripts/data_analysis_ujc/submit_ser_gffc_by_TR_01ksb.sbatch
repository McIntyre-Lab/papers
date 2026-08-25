#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gffcK
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/gffcK_%A-%a.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --account=mcintyre
##SBATCH --qos=mcintyre-b
#SBATCH --array=2-49

# gff compare on TRs

module purge
module load gffcompare/0.12.2

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
ALN=${PROJ}/aln_ser_dorReads_2_dser1

UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

OUT=$PROJ/gffcompare_by_TR
	mkdir -p $OUT

## Get design file information (2 - 49)
DESIGN_FILE=$PROJ/design_files/design_file_sample2BC_w_origDataPath_04amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
BC=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[3]}
TECH=${ARRAY[4]}
DATE=${ARRAY[5]}
ONT_DIR=${ARRAY[6]}

SAMPLEID=${SAMPLE}_TR${TR}

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dser1_ujc.gtf

echo "Assigning reference genes..."

GTF=${ALN}/${SAMPLEID}_2_dser1_noMQ0_uniq.gtf

# setting prefix for gffcompare output
NAME=${SAMPLEID}_2_dser1_noMQ0_uniq

date

## run GFFcompare
gffcompare -r $REF_GTF -o $OUT/${NAME}_roz $ALN/${NAME}.gtf

## use avn script to  update starting gtf with geneIDs
python $UTIL/correct_gffcompare_GTF_gene_id.py \
-a $OUT/${NAME}_roz.annotated.gtf \
-g $ALN/${NAME}.gtf \
-o $OUT/${SAMPLEID}_2_dser1_updGeneID.gtf \
--key ${OUT}/${SAMPLEID}_2_dser1_gffcompare_gene_key.csv \
--bin-TS ${OUT}/${SAMPLEID}_2_dser1_potential_trans_spliced_reads.csv

date

rm $OUT/$NAME*roz*
rm $ALN/$NAME*tmap
rm $ALN/$NAME*refmap
