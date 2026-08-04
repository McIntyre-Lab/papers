#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sqK
#SBATCH --qos=mcintyre
#SBATCH --account=mcintyre
#SBATCH --output=/blue/mcintyre/share/kopp_lmm_head_data/scripts/SLURM_LOGS/sqK_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=3-49

## run sqanti with UJC gtf for dmel650


## input dirs
PROJ=/blue/mcintyre/share/kopp_lmm_head_data
PYCHOP=${PROJ}/dorado_basecalling/pychopper
ALN=${PROJ}/aln_ser_dorReads_2_dser1

SQANTI=${PROJ}/sqanti3QC_ser_dorReads_dser1_ujc
    mkdir -p $SQANTI

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

ROZ=$PROJ/roz_sq_${SAMPLEID}
    mkdir -p $ROZ


    SAM=${ALN}/${SAMPLEID}_fl_rsc_2_dser1.sam
    BAM=${ALN}/${SAMPLEID}_fl_rsc_2_dser1.bam
    GTF=${ALN}/${SAMPLEID}_fl_rsc_2_dser1.gtf
    BED=${ALN}/${SAMPLEID}_fl_rsc_2_dser1.bed

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

    REF_GTF=/blue/mcintyre/share/sex_specific_splicing/within_genome_uniq_ujc/dser1_all_uniq_ujc_corrected_associated_gene.gtf
    REF_FA=/blue/mcintyre/share/references/dser1.1/GCF_002093755.2/GCF_002093755.2_Dser1.1_genomic.fna

    module purge
    module load sqanti3/4.2

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID}_fl_rsc_dser1 \
        ${GTF} \
        ${REF_GTF} \
        ${REF_FA}

    ## Move all files to output directory
    mv ${ROZ}/*.* ${SQANTI}/.

## Get classification counts

    MAPPED=$(awk 'NR!=1' ${SQANTI}/${SAMPLEID}_fl_rsc_dser1_classification.txt | wc -l)

    awk -F "\t" 'NR!=1{print $6}' ${SQANTI}/${SAMPLEID}_fl_rsc_dser1_classification.txt | \
        sort | uniq -c | \
         awk -v total=${MAPPED} 'BEGIN{print "structural_category,count,proportion"}{print $2","$1","$1/total}' \
        > ${SQANTI}/${SAMPLEID}_fl_rsc_dser1_class_cnt.csv

done

rm -r $ROZ

