#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=upd_genome
#SBATCH --output=/blue/mcintyre/share/pxs_lmm_dros_data/scripts/update_genome/SLURM_LOGS/upd_genome_%A_%a.out
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --array=1,2,4

## update w1118 for each RIL - testing with 12272


module load bcftools/1.19 bwa/0.7.17 samtools/1.18 python/3.6

## Set directories
PROJ=/blue/mcintyre/share/pxs_lmm_dros_data
IND=$PROJ/aln_upd_genome_fl_rsc_unclass
OUTD=$PROJ/updated_w1118_tester_genomes
    mkdir -p $OUTD

## Get design file information
DESIGN_FILE=$PROJ/design_files/list_rils.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RIL=${ARRAY[0]}

echo "${RIL}
"
date

## Get reference FASTA file to update 
FASTA=/blue/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
#FASTA=$PROJ/updated_genomes/w1118_snp_upd_genome_chromSubset.fasta

#VCF=${IND}/W1118_snps_w_ril_12272_snps.vcf.gz (manually added additional header info)
VCF=${IND}/W1118_snp_upd_${RIL}_snp.vcf.gz
NAME=w1118_${RIL}_snps

## Update genome
echo "Updating genome...
"
cat ${FASTA} | bcftools consensus ${VCF} \
    > ${OUTD}/${NAME}.fasta

date

## Build BWA and samtools indeces for updated genome
echo "build bwa index
"
bwa index -p ${OUTD}/${NAME}_BWA -a bwtsw ${OUTD}/${NAME}.fasta
echo "build samtools faidx
"
samtools faidx ${OUTD}/${NAME}.fasta


## Check updated genome
echo "Check updated genome...
"
## Convert both reference and updated fa files to single line fa (seq on 1 line)
if [[ ! -e $(dirname ${FASTA})/$(basename ${FASTA} .fasta).oneLine.fasta ]]; then
    awk '/^>/{printf("\n%s\n",$0);next;}{printf("%s",$0);} \
        END{printf("\n");}' ${FASTA} | tail -n +2 \
        > $(dirname ${FASTA})/$(basename ${FASTA} .fasta).oneLine.fasta
fi

awk '/^>/{printf("\n%s\n",$0);next;}{printf("%s",$0);} \
    END{printf("\n");}' ${OUTD}/${NAME}.fasta | tail -n +2 \
    > ${OUTD}/${NAME}.oneLine.fasta

diff ${OUTD}/${NAME}.oneLine.fasta \
    $(dirname ${FASTA})/$(basename ${FASTA} .fasta).oneLine.fasta \
    > ${OUTD}/${NAME}_check_diff.txt


rm ${OUTD}/${NAME}.oneLine.fasta

## Build BWA and samtools indeces for updated genome chromsome subset
echo "build bwa index
"
bwa index -p ${OUTD}/${NAME}_BWA -a bwtsw ${OUTD}/${NAME}.fasta
echo "build samtools faidx
"
samtools faidx ${OUTD}/${NAME}.fasta

