#!/bin/bash
#SBATCH --mail-user=k.bankole@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gffcT
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/gffcT_%A-%a.out
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45gb
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --array=2-9

# gff compare on TRs

module purge
module load gffcompare/0.12.2

export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
UTIL=/blue/mcintyre/k.bankole/github/TranD/utilities

## Get design file information (2 - 9)
DESIGN_FILE=$PROJ/design_files/design_dyak_dsan_parent_2_dir_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FLAG=${ARRAY[0]}
SAMPLE=${ARRAY[1]}
BC=${ARRAY[2]}
BCNUM=${ARRAY[3]}
TR=${ARRAY[4]}
TECH=${ARRAY[5]}
DATE=${ARRAY[6]}
ONT_DIR=${ARRAY[7]}

SAMPLEID=${SAMPLE}_TR${TR}

ALN=${PROJ}/aln_dorReads_2_dyak2_or_dsan1

if [[ $SAMPLEID == dsan* ]] ; then
        GENOME=dsan1
else
    	GENOME=dyak2
fi

OUT=$PROJ/gffcompare_by_TR
	mkdir -p $OUT

REF_GTF=/blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_${GENOME}_ujc.gtf

echo "Assigning reference genes..."

GTF=${ALN}/${SAMPLEID}_2_${GENOME}_noMQ0_uniq.gtf

# setting prefix for gffcompare output
NAME=${SAMPLEID}_2_${GENOME}_noMQ0_uniq

date

## run GFFcompare
gffcompare -r $REF_GTF -o $OUT/${NAME}_roz $ALN/${NAME}.gtf

## use avn script to  update starting gtf with geneIDs
python $UTIL/correct_gffcompare_GTF_gene_id.py \
-a $OUT/${NAME}_roz.annotated.gtf \
-g $ALN/${NAME}.gtf \
-o $OUT/${SAMPLEID}_2_${GENOME}_updGeneID.gtf \
--key ${OUT}/${SAMPLEID}_2_${GENOME}_gffcompare_gene_key.csv \
--bin-TS ${OUT}/${SAMPLEID}_2_${GENOME}_potential_trans_spliced_reads.csv

date

rm $OUT/$NAME*roz*
rm $ALN/$NAME*tmap
rm $ALN/$NAME*refmap


