#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/sqanti_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-11

### SQANTI QC of maize ignored clusters following tofu2 – using passing clusters + B73 ref as SQANTI reference

#module load pacbio/5.1.0 isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0
## pacbio/5.1.0 module no longer available, may not need
module load isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0
module load bedtools/2.29.2

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/compare_b73_2_mo17/sqanti_b73_ignored_clusters
    mkdir -p ${OUTD}
TOFU=$PROJ/compare_b73_2_mo17/tofu2_b73
MAP=$PROJ/compare_b73_2_mo17/mapping_minimap2_b73
IGNORE=$PROJ/compare_b73_2_mo17/uncollapsed_b73_IsoSeq_reads

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

REF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
REF_GTF=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf

## Set input and output directories
OUT=${OUTD}/${SAMPLE_NAME}
    mkdir -p ${OUT}
INTOFU=${TOFU}/${SAMPLE_NAME}
INMAP=${MAP}/${GENO}/${TRT}
ROZ=${OUT}/roz_sqanti
    mkdir -p ${ROZ}

## Loop over ignored types
for TYPE in low_cov low_ident; do

    ## Subset SAM file for ignored cluster IDs
    python $PROJ/scripts/Common_extracFeatureFromList.py \
        -l ${IGNORE}/list_${SAMPLE_NAME}.ignored_${TYPE}_ids.b73_ref.txt \
        -i ${INMAP}/${SAMPLE_NAME}.polished.all.hq.mapped.sam \
        -t sam \
        -o ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.sam

    ## Convert sam to bam
    samtools view -b ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.sam \
        > ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.bam

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
        -i ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.bam \
        > ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.bed
    awk -v source=${SOURCE} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' \
        ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.bed \
        > ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.gtf

    ## GTF (is GTF but has .gff ending) from output of collapse_isoforms_by_sam.py
    PACBIO_GTF=${INTOFU}/${SAMPLE_NAME}.collapsed.gff

    ## Concatenate the Reference GTF and collapsed cluster GTF
    cat ${PACBIO_GTF} ${REF_GTF} \
        > ${ROZ}/collapse_cluster_plue_ref.gtf

    date
    echo "***SQANTI***
"
    python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py \
        -g \
        ${ROZ}/${SAMPLE_NAME}.ignored_${TYPE}.gtf \
        ${ROZ}/collapse_cluster_plue_ref.gtf \
        ${REF} -o ${SAMPLE_NAME}_${TYPE} -d ${OUT} \
        2>${OUT}/${SAMPLE_NAME}_${TYPE}.sqanti.log

    ## Remove uneccesary output to limit space
    rm -r ${OUT}/GMST
    rm -r ${OUT}/RTS
    rm ${OUT}/${SAMPLE_NAME}.ignored_${TYPE}_corrected.faa
    rm ${OUT}/${SAMPLE_NAME}.ignored_${TYPE}_corrected.fasta

done

#rm -r ${ROZ}
