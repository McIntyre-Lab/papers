#!/bin/bash
#SBATCH --mail-user=ashley.myrick@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH--job-name=Dorsq
#SBATCH --qos=mcintyre-b
#SBATCH --account=mcintyre
#SBATCH--output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/Dorsq_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --array=2-68
##SBATHC --array=2

### running QC  on 01h_19july2022 run
## sqanti3 QC on  1) fl + res and 2) unclass

PROJ=/blue/mcintyre/share/pxs_lmm_dros_data

PYCHOP=${PROJ}/dorado_basecalling/pychopper
ALN=${PROJ}/dorado_basecalling/aln_genome_minimap2_dmel6

SQANTI=${PROJ}/dorado_basecalling/sqanti3QC_ref_dmel650_dsimWXD1_union
    mkdir -p $SQANTI

REF_FA=/blue/mcintyre/share/references/dmel_fb650/dmel-all-chromosome-r6.50.fasta
REF_GTF=/blue/mcintyre/share/sex_specific_splicing/dm650_vs_dsWXD1_TMM_dm650/dm650_2_dm650_dmWXD_2_dm650_union_UJC.gtf

## Get design file information (2 - 109)
DESIGN_FILE=$PROJ/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
SAMPLEID=${ARRAY[1]}
BCNUM=${ARRAY[2]}
TR=${ARRAY[8]}
TECH=${ARRAY[9]}
DATE=${ARRAY[10]}
BCNUM=${ARRAY[13]}

ROZ=${PROJ}/roz_sq_${SAMPLEID}
    mkdir -p ${ROZ}

module load samtools/1.10
module load bedtools/2.29.2

for TYPE in fl_rsc unclass ; do

    MAPPED=${ALN}/${SAMPLEID}_${TYPE}_minimap2.sam
    BAM=${ALN}/${SAMPLEID}_${TYPE}_minimap2.bam
    GTF=${ALN}/${SAMPLEID}_${TYPE}_minimap2.gtf
    BED=${ALN}/${SAMPLEID}_${TYPE}_minimap2.bed

    ## SAM to BAM
    samtools view -b -F 2048 ${MAPPED} > ${BAM}

    ## BAM to BED
    bedtools bamtobed -split -i ${BAM} > ${BED}

    ## BED to GTF
    SOURCE=${SAMPLEID}_${TYPE}
    awk -v source=${SOURCE} -F "\t" '{print $1"\t"source"\texon\t"$2+1"\t"$3"\t.\t"$6"\t.\tgene_id \""$4"\"; transcript_id \""$4"\";"}' $BED > $GTF

    ## Remove intermediate files
    rm ${BED} ${BAM}


    module purge
    module load sqanti3/4.2

        ##  NO --gtf option (Now is default and not required)
        ##  --skipORF option to skip ORF prediction (runs much faster)
        ##  --report pdf option to make a SQANTI QC report in PDF format (HTML fails in current version)
        ##  --force_id_ignore option to allow transcript ID names that are not like PB.X.Y
        ##  Order of files is (1) isoforms (FA/FQ/GTF), (2) annotation (ref GTF/GFF), (3) genome (ref FA)

        #   --dir is path to output directory
        #   -o prefix name to use

        ## order of files input into sqanti3_qc.py is critical!!
        # 1. sam2gtf (from read align to genome)
        # 2. FSM consolidated ref gtf file 
        # 3. ref genome fasta file

        ## Sqanti3 QC

    sqanti3_qc.py \
        --skipORF \
        --report pdf \
        --force_id_ignore \
        --dir ${ROZ} \
        -o ${SAMPLEID}_${TYPE} \
        ${ALN}/${SAMPLEID}_${TYPE}_minimap2.gtf \
        ${REF_GTF} \
        ${REF_FA}

    ## Move all files to output directory
    mv ${ROZ}/*.* ${SQANTI}/.


## Get classification counts

    if [[ ${TYPE} == "fl_rsc" ]]; then
        INFQ=${ROZ}/${SAMPLEID}_${TYPE}.fq
        cat \
            ${PYCHOP}/fl_${SAMPLEID}.fq \
            ${PYCHOP}/rescue_${SAMPLEID}.fq 2>/dev/null > ${INFQ}
    elif [[ ${TYPE} == "unclass" ]]; then
        INFQ=${PYCHOP}/unclass_${SAMPLEID}.fq
    else
        echo "can't find fq for mapped counting"
    fi

    TOTAL=$(awk 'END {print NR/4}' ${INFQ})
    MAPPED=$(awk 'NR!=1' ${SQANTI}/${SAMPLEID}_${TYPE}_classification.txt | wc -l)

    awk -F "\t" 'NR!=1{print $6}' ${SQANTI}/${SAMPLEID}_${TYPE}_classification.txt | \
        sort | uniq -c | \
         awk -v total=${MAPPED} 'BEGIN{print "structural_category,count,proportion"}{print $2","$1","$1/total}' \
        > ${SQANTI}/${SAMPLEID}_${TYPE}_class_cnt.csv

    ## Remove unnecessary files
    rm ${SQANTI}/${SAMPLEID}_${TYPE}_corrected.fasta
    rm ${SQANTI}/${SAMPLEID}_${TYPE}_corrected.genePred
    rm ${SQANTI}/${SAMPLEID}_${TYPE}_corrected.gtf
    rm ${SQANTI}/${SAMPLEID}_${TYPE}_corrected.gtf.cds.gff
    rm ${SQANTI}/refAnnotation_${SAMPLEID}_${TYPE}.genePred

done

rm -r $ROZ
