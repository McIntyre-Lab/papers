#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=aln_upd
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/SLURM_LOGS/aln_upd_%A_%a.out
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --array=1,2,4

## align all reps to upd RIL

module load minimap/2.12
module load samtools/1.12
module load gcc/5.2.0
module load bioawk/1.0
module load bcftools
module load bedtools 

### Set Project Directory
PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
PYCHOP=$PROJ/dorado_basecalling/pychopper
VCF=$PROJ/snp_tables_by_genotype

OUTD=$PROJ/aln_upd_genome_fl_rsc_unclass
    mkdir -p $OUTD

## Get design file information ( 1-4 )
DESIGN_FILE=$PROJ/design_files/list_rils.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RIL=${ARRAY[0]}

ROZ=${PROJ}/roz_${RIL}
    mkdir -p ${ROZ}

## Concatenate full length (fl) and rescue (rsc) and unclass FASTQ files for the given RIL
cat ${PYCHOP}/*${RIL}*.fq > ${ROZ}/dm${RIL}.fastq

REF=$PROJ/updated_genomes_0.5_rmRefMismatch/ril_${RIL}_snp_0.5_upd_genome.fasta
REFNAME=upd_${RIL}

## Do not use the -u f option with unclassified reads since
##     they have not had the orientation correction like fl+rs
minimap2 \
    -t 8 \
    -a \
    -x splice \
    --secondary=no \
    -C 5 \
    ${REF} \
    ${ROZ}/dm${RIL}.fastq \
    > ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome.sam

## remove reads with mapping quality of 0
samtools view \
   -h \
   -q 1 \
    -bo \
    ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_noMQ0.bam \
    ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome.sam

## include reads without flag 2052 (omit unmapped and supplementary)
samtools view \
    -h \
    -F 2052 \
    ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_noMQ0.bam \
    > ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_uniq.bam

# sort
samtools sort -o ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_uniq.sorted.bam \
    ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_uniq.bam

# create vcf file
bcftools mpileup -f $REF ${OUTD}/dm${RIL}_2_${REFNAME}_upd_genome_uniq.sorted.bam | \
bcftools call -m -o $OUTD/dm${RIL}_2_${REFNAME}_upd_genome.vcf

# filter vcf - keep with more than 5 reads
bcftools filter -i  'DP>5' $OUTD/dm${RIL}_2_${REFNAME}_upd_genome.vcf -o $OUTD/dm${RIL}_2_${REFNAME}_upd_genome_DP5.vcf

## compare with bedtools - above vcf and starting vcf - output variants in both
bedtools intersect -a $OUTD/dm${RIL}_2_${REFNAME}_upd_genome_DP5.vcf \
    -b $VCF/ril_${RIL}_snp_filter_0.5_rmRefMismatch.vcf \
    > $OUTD/dm${RIL}_2_${REFNAME}_int_noH_0.5_rmRefMismatch.vcf

## keep snps 
awk -F'\t' 'BEGIN {OFS="\t"} {if (length($4) == 1 && length($5) == 1 || /^#/) print}' $OUTD/dm${RIL}_2_${REFNAME}_int_noH_0.5_rmRefMismatch.vcf \
    > $OUTD/dm${RIL}_2_${REFNAME}_int_noH_0.5_rmRefMismatch_snps.vcf

# add header
sed "s/CHANGEME/dm${RIL}/g" $OUTD/header.txt > $OUTD/header_${RIL}.txt
cat $OUTD/header_${RIL}.txt $OUTD/dm${RIL}_2_${REFNAME}_int_noH_0.5_rmRefMismatch_snps.vcf > $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf

## create compressed vcf and index
bgzip -f $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf > $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf.gz
bcftools index -f $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf.gz -o $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf.gz.csi

## only GT=0/0
bcftools view -i 'GT="0/0"' $OUTD/dm${RIL}_2_${REFNAME}_int_0.5_rmRefMismatch_snps.vcf.gz -o $OUTD/dm${RIL}_2_${REFNAME}_upd_genome_gt00.vcf.gz

## pull these postions from starting RIL vcf file
bedtools intersect -a $VCF/ril_${RIL}_snp_filter_0.5_rmRefMismatch.vcf \
    -b $OUTD/dm${RIL}_2_${REFNAME}_upd_genome_gt00.vcf.gz \
    > $OUTD/dm${RIL}_snps_gt00_w1118.vcf

## cat to w1118 vcf file for updating
cat $VCF/W1118_snp_renamed_sorted.vcf $OUTD/dm${RIL}_snps_gt00_w1118.vcf > $OUTD/W1118_snp_w_ril_${RIL}.vcf

## sort and remove duplicate rows based on chrom position ref and alt columns
grep -v '^#' $OUTD/W1118_snp_w_ril_${RIL}.vcf | sort -k1,1 -k2,2n -k4,4 -k5,5 | awk '!seen[$1, $2, $4, $5]++' >> $OUTD/W1118_snp_w_ril_${RIL}_sort.vcf

sort -k1,1 -k2,2n $OUTD/W1118_snp_w_ril_${RIL}_sort.vcf > $OUTD/W1118_snp_w_ril_${RIL}_sort2.vcf

grep '^#' $OUTD/W1118_snp_w_ril_${RIL}.vcf > $OUTD/W1118_snp_w_ril_${RIL}_header.vcf

cat $OUTD/W1118_snp_w_ril_${RIL}_header.vcf \
    $OUTD/W1118_snp_w_ril_${RIL}_sort2.vcf \
    > $OUTD/W1118_snp_upd_${RIL}_snp.vcf

bgzip -f $OUTD/W1118_snp_upd_${RIL}_snp.vcf
tabix -f $OUTD/W1118_snp_upd_${RIL}_snp.vcf.gz


