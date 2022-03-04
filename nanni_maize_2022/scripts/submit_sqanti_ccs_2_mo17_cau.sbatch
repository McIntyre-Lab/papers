#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/sqanti_ccs_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-11

### SQANTI QC of maize mapped ccs reads – using Mo17 CAU ref as SQANTI reference

module load isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0 bedtools/2.29.2

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/compare_b73_2_mo17/sqanti_mo17_cau_ccs_reads
    mkdir -p ${OUTD}
IND=$PROJ/check_isoseq3_lost_genes/minimap2

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

REF=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa
REF_GTF=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0_mgdb.gtf

## Set input and output directories
OUT=${OUTD}/${SAMPLE_NAME}
    mkdir -p ${OUT}
ROZ=${OUT}/roz_sqanti
    mkdir -p ${ROZ}

## Convert SAM to GTF using SAM to BED12 –> BED12 to GTF

## Convert sam to bam
samtools view -b ${IND}/${SAMPLE_NAME}_minimap2_mo17_cau.sam \
    > ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.bam

## Report each part of split bam alignment as distinct BED12 interval
## Create GTF file with the following columns:
##     chr, source, feature, start, end, score, strand, frame, attribute
## source -  unique label indicating where the annotations came from
## feature - will only be exon in this conversion
## start - 1-based
## score - degree of confidence in the feature's existence and coordinate, not necessary so set as "."
## frame - indicator of feature codon frame, not necessary so set as "."
## attribute - includes the two mandatory attributes "gene_id" and "transcript_id"
bedtools bamtobed \
    -split \
    -i ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.bam \
    > ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.bed
awk -v source=${SAMPLE_NAME} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' \
    ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.bed \
    > ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.gtf

date
echo "***SQANTI***
"
python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py \
    -g \
    ${ROZ}/${SAMPLE_NAME}_minimap2_mo17_cau.gtf \
    ${REF_GTF} \
    ${REF} -o ${SAMPLE_NAME}_ccs -d ${OUT} \
    2>${OUT}/${SAMPLE_NAME}_ccs.sqanti.log

## Remove uneccesary output to limit space
rm -r ${OUT}/GMST
rm -r ${OUT}/RTS
rm ${OUT}/${SAMPLE_NAME}_minimap2_mo17_cau_corrected.faa
rm ${OUT}/${SAMPLE_NAME}_minimap2_mo17_cau_corrected.fasta

rm -r ${ROZ}
