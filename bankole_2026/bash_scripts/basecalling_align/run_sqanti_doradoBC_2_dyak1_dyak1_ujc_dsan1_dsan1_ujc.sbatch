#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sqTO
#SBATCH --qos=mcintyre-cid
#SBATCH --account=mcintyre-cid
#SBATCH --output=/blue/mcintyre/share/transcript_ortholog/scripts/SLURM_LOGS/sqTO_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=86gb
##SBATCH --array=3
#SBATCH --array=9

## sqanti3QC and class count files

## input dirs
PROJ=/blue/mcintyre/share/transcript_ortholog
ALN=${PROJ}/aln_dorReads_2_dyak1_or_dsan1

SQANTI=${PROJ}/sqanti3QC_yak_san_dorReads_dyak1_dsan1_ujc
    mkdir -p $SQANTI

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

ROZ=$PROJ/roz_sq_${SAMPLEID}
    mkdir -p $ROZ


if [[ $SAMPLEID == dyak* ]] ; then
    SAM=${ALN}/${SAMPLEID}_fl_rsc_2_dyak1.sam
    BED=${ALN}/${SAMPLEID}_fl_rsc_2_dyak1.bed
    GTF=${ALN}/${SAMPLEID}_fl_rsc_2_dyak1.gtf
    BAM=${ALN}/${SAMPLEID}_fl_rsc_2_dyak1.bam
    BAMS=${ALN}/${SAMPLEID}_fl_rsc_2_dyak1_sort.bam
    TAG=dyak1
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/within_genome_uniq_ujc/dyak1_all_uniq_ujc_corrected_associated_gene.gtf
    REF_FA=/blue/mcintyre/share/references/dyak_rogers/Dyak-all-chrom-wb-r1.05.fa
else
    SAM=${ALN}/${SAMPLEID}_fl_rsc_2_dsan1.sam
    BED=${ALN}/${SAMPLEID}_fl_rsc_2_dsan1.bed
    GTF=${ALN}/${SAMPLEID}_fl_rsc_2_dsan1.gtf
    BAM=${ALN}/${SAMPLEID}_fl_rsc_2_dsan1.bam
    BAMS=${ALN}/${SAMPLEID}_fl_rsc_2_dsan1_sort.bsam
    TAG=dsan1
    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/within_genome_uniq_ujc/dsan1_all_uniq_ujc_corrected_associated_gene.gtf
    REF_FA=/blue/mcintyre/share/references/dsan_Prin_1.1/GCF_016746245.2/GCF_016746245.2_Prin_Dsan_1.1_genomic.fna
fi


module purge
module load samtools/1.10 bedtools/2.29.2

## SAM to BAM

samtools view -b -F 2048 ${SAM} > ${BAM}

## BAM to BED
bedtools bamtobed -split -i ${BAM} > ${BED}

## BED to GTF
SOURCE=${SAMPLEID}_fl_rsc
awk -v source=${SOURCE} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

## Remove intermediate files
rm ${BED} ${BAM}


    module purge
    module load sqanti3/4.2

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID}_fl_rsc_${TAG} \
        ${GTF} \
        ${REF_GTF} \
        ${REF_FA}

    ## Move all files to output directory
    mv ${ROZ}/*.* ${SQANTI}/.


## Get classification counts

    MAPPED=$(awk 'NR!=1' ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt | wc -l)

    awk -F "\t" 'NR!=1{print $6}' ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_classification.txt | \
        sort | uniq -c | \
         awk -v total=${MAPPED} 'BEGIN{print "structural_category,count,proportion"}{print $2","$1","$1/total}' \
        > ${SQANTI}/${SAMPLEID}_fl_rsc_${TAG}_class_cnt.csv

done

rm -r $ROZ

