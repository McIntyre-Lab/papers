#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=aln_upd
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/aln_upd_%A_%a.out
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=26gb
#SBATCH --array=2-68
##SBATCH --array=2

module load minimap/2.12
module load samtools/1.12
module load python/2.7.6
module load bedtools/2.30.0

### Set Project Directory
PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
PYCHOP=$PROJ/dorado_basecalling/pychopper
OUTD=$PROJ/aln_upd_gene_02amm
    mkdir -p $OUTD
SAMC=$PROJ/ase_counts_2_upd_w1118_and_2_upd_ril_02amm
    mkdir -p $SAMC


DESIGN_FILE=$PROJ/design_files/sample_sampleID_bc_dsgn_w_origDataPath_02amm.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE=${ARRAY[0]}
SAMPLEID=${ARRAY[1]}
GENO=${ARRAY[3]}
TRT=${ARRAY[4]}
TR=${ARRAY[8]}
TECH=${ARRAY[9]}
DATE=${ARRAY[10]}
BCNUM=${ARRAY[13]}

RIL=${GENO:2}

ROZ=${PROJ}/roz_${SAMPLEID}
    mkdir -p ${ROZ}

## Concatenate full length (fl) and rescue (rsc) and unclass FASTQ files for the given sample
##     and fix FQ names (remove all extra info from ID line other than name)
##     samcompare is not able to split the other info out
cat ${PYCHOP}/fl_${SAMPLEID}.fq \
    ${PYCHOP}/rescue_${SAMPLEID}.fq \
    ${PYCHOP}/unclass_${SAMPLEID}.fq \
    > ${ROZ}/${SAMPLEID}.fastq

## Fix fastq name to be shorter for samcompare
awk '{if(NR%4==1){print $1}else{print}}' ${ROZ}/${SAMPLEID}.fastq \
    > ${ROZ}/${SAMPLEID}_fixName.fastq

## Get average read length in FQ file
AVE_RL=$(awk '{if(NR%4==2 && length) {count++; bases += length} } END {print bases/count}' ${ROZ}/${SAMPLEID}.fastq | awk '{ printf "%.0f",$1 }')
echo "Average RL = ${AVE_RL}"

# get gene FA for updated RIL
URIL=$PROJ/updated_genomes_0.5_rmRefMismatch/ril_${RIL}_snp_0.5_upd_genome_gene_plusMinus100.fasta

# get gene FA for new  updated TESTER (w1118)
UTEST=$PROJ/updated_w1118_tester_genomes/w1118_${RIL}_snp_upd_genome_gene_plusMinus100.fasta

## align to both updated w1118  and to updated ril
for WHICH in line tester; do
    if [[ ${WHICH} == "line" ]]; then
        REF=${URIL}
        REFNAME=${GENO}
    else
        REF=${UTEST}
        REFNAME=w1118
   fi

    minimap2 \
        -t 16 \
        -a \
        -x splice \
        --secondary=yes \
        -N 200 \
        -p 0.9 \
        -C 5 \
        ${REF} \
        ${ROZ}/${SAMPLEID}.fastq \
        > ${OUTD}/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100.sam

    ## remove reads with mapping quality of 0
    samtools view \
        -h \
        -q 1 \
        ${OUTD}/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100.sam \
        > ${OUTD}/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100_noMQ0.sam

    ## exclude header, include reads without flag 2052 (omit unmapped and supplementary)
    samtools view \
        -F 2052 \
        ${OUTD}/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100_noMQ0.sam \
        > ${OUTD}/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100_noMQ0_uniq.sam

    #rm ${OUTD}/${SAMPLE_NAME}_2_${REFNAME}_upd_gene_plusMinus100.sam

    ## count reads to gene
    (echo "geneID,count"; awk '{print $3}' $OUTD/${SAMPLEID}_2_${REFNAME}_upd_gene_plusMinus100_noMQ0_uniq.sam | sort | uniq -c | sort -nr | awk '{print $2","$1}') \
        > $OUTD/${SAMPLE}_2_${REFNAME}_gene_cnts.csv

#done

## run samcompare

## Get sam files
SAMRIL=${OUTD}/${SAMPLEID}_2_${GENO}_upd_gene_plusMinus100_noMQ0_uniq.sam
SAMTEST=${OUTD}/${SAMPLEID}_2_w1118_upd_gene_plusMinus100_noMQ0_uniq.sam

## Reorder BED file
GBED=/blue/mcintyre/share/references/dmel_fb643/dmel643_gene_region_plus_minus_100.bed
awk -v OFS='\t' '{print $4,$2,$3,$1}' ${GBED} \
    > ${ROZ}/gene_bed_reorder.bed

# create direcory for reads id output of samcompare categories
mkdir -p ${SAMC}/readIDs

python $PROJ/scripts/sam_compare_w_feature_reads_out.py  \
    -l $AVE_RL \
    -f ${ROZ}/gene_bed_reorder.bed \
    -q ${ROZ}/${SAMPLEID}_fixName.fastq \
    -B ${SAMTEST} \
    -A ${SAMRIL} \
    -g ${SAMC}/samcompare_${SAMPLEID}.log \
    -c ${SAMC}/ase_counts_${SAMPLEID}.csv \
    -t ${SAMC}/ase_totals_${SAMPLE}.txt \
    -n \
    -p ${SAMC}/readIDs_

done 

rm -r ${ROZ}

