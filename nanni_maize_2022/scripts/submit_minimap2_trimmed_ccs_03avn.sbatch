#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=minimap2
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/minimap2_%A_%a.out
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11

## Map to maize B73 genome with minimap2

module load minimap/2.12 samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/check_isoseq3_lost_genes/minimap2
    mkdir -p ${OUTD}
LOG=${OUTD}/log
    mkdir -p ${LOG}
IND=$PROJ/check_isoseq3_lost_genes/bam2fq

## Get reference genome FASTA
B73=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
CAU=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa
YAN=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.fasta

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

date

for NAME in b73 mo17_cau mo17_yan; do
    if [[ ${NAME} == "b73" ]]; then
        REF=${B73}
    elif [[ ${NAME} == "mo17_cau" ]]; then
        REF=${CAU}
    elif [[ ${NAME} == "mo17_yan" ]]; then
        REF=${YAN}
    else
        echo "OOPS"
    fi
    echo "

*** Align to ${GENO} genome with minimap2 ***

Input: ${IND}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.fq
REF:   ${REF}
"
#    minimap2 \
#        -t 12 \
#        -a \
#        -x splice \
#        -u f \
#        --secondary=no \
#        -C 5 \
#        ${REF} \
#        ${IND}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.fq \
#        > ${OUTD}/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.sam 2>$LOG/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.log

    ## include header, include reads without flag 2052 (omit unmapped and supplementary)
#    samtools view \
#        -h \
#        -F 2052 \
#        ${OUTD}/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.sam \
#        > ${OUTD}/${SAMPLE_NAME}_minimap2_${NAME}.sam

    date

    ## Counting mapped reads, unmapped reads, and multi-mapped reads (has supp alignments - sam flag 2048)
    echo "sampleID,ref,num_start_read,num_mapped_read,num_unmapped_read,sum_map_unmap_read,num_multimap_read,prop_mapped_read,prop_unmap_read,prop_multimap_read" \
        > ${OUTD}/${SAMPLE_NAME}_2_${NAME}_ccs_mapped_read_count.csv

    MAP=$(samtools view ${OUTD}/${SAMPLE_NAME}_minimap2_${NAME}.sam | wc -l)

    UNMAP=$(samtools view -f 4 ${OUTD}/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.sam | wc -l)

    MMAP=$(samtools view -f 2048 ${OUTD}/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.sam | cut -f 1 | sort | uniq | wc -l)

    awk -v ref=${NAME} -v map=${MAP} -v unmap=${UNMAP} -v multi=${MMAP} -v sampleid=${SAMPLE_NAME} -v OFS="," \
        'END{total=NR/4; print sampleid, ref, total, map, unmap, map+unmap, multi, map/total, unmap/total, multi/total}' \
        ${IND}/${SAMPLE_NAME}.ccs.trimmed.primer_5p--primer_3p.fq \
        >> ${OUTD}/${SAMPLE_NAME}_2_${NAME}_ccs_mapped_read_count.csv

#    rm ${OUTD}/${SAMPLE_NAME}_unfiltered_minimap2_${NAME}.sam

done
