#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/sqanti_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-11

### SQANTI isoseq QC of maize PacBio data -- Mo17 CAU ref

module load pacbio/5.1.0 isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/sqanti_mo17_cau
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IND=$PROJ/tofu2_mo17_cau

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}
REF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa

REF_GTF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0_mgdb.gtf

## Set input and output directories to the correct species
OUT=$OUTD/$GENO/$TRT
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IN=$IND/$GENO/$TRT

## GTF (is GTF but has .gff ending) from output of collapse_isoforms_by_sam.py
PACBIO_GTF=$IN/${SAMPLE_NAME}.collapsed.gff

ABUNDANCE=$IN/${SAMPLE_NAME}.collapsed.abundance.txt

date
echo "***SQANTI***
"
python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py -g ${PACBIO_GTF} ${REF_GTF} \
    $REF -o ${SAMPLE_NAME} -d $OUT -fl $ABUNDANCE 2>$OUT/${SAMPLE_NAME}.sqanti.log

