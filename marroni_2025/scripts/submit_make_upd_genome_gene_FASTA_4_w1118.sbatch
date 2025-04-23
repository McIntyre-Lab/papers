#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=geneFA
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/update_genome/SLURM_LOGS/gene_upd_FA_%A.out
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --array=1,2,4

module load bedtools/2.30.0
module load samtools/1.15
module load bwa/0.7.17

## Make gene FASTA for updated w1118 genomes

### Set Directories
PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
IND=/blue/mcintyre/share/references/dmel_fb643
TESTER=$PROJ/updated_w1118_tester_genomes

## Get dmel 6.43 gene BED file (+/-100nt)
BED=${IND}/dmel643_gene_region_plus_minus_100.bed

## Get design file information
DESIGN_FILE=$PROJ/design_files/list_rils.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RIL=${ARRAY[0]}

echo "${RIL}
"
date

## Get gene FASTA for w1118
## Used -name argument to use gene_id for FASTA name field
bedtools getfasta \
    -fi ${TESTER}/w1118_${RIL}_snps.fasta \
    -bed ${BED} \
    -nameOnly \
    -fo ${TESTER}/w1118_${RIL}_snp_upd_genome_gene_plusMinus100.fasta

samtools faidx ${TESTER}/w1118_${RIL}_snp_upd_genome_gene_plusMinus100.fasta

bwa index \
    -p ${TESTER}/w1118_${RIL}_snp_upd_genome_gene_plusMinus100_BWA \
    -a bwtsw \
    ${TESTER}/w1118_${RIL}_snp_upd_genome_gene_plusMinus100.fasta

