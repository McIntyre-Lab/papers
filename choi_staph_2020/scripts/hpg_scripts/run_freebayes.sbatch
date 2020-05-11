#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --job-name=callSnps
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/staph_relapse/scripts/SLURM_LOGS/callSnps_%A-%a.out
#SBATCH --time=6:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb    		# per processor memory
#SBATCH --array=3-79


#### call snps with freebayes
## HAPLOID!!!
##	-p 1

module load samtools/0.1.19 python freebayes/0.9.15 vcflib vcftools

## Set directories
PROJ=/ufrc/mcintyre/share/staph_relapse
REF_DIR=${PROJ}/CC_references
ALN=${PROJ}/aln_bwa_PE

## Create output directory
OUTPUT=${PROJ}/freebayes
    mkdir -p ${OUTPUT}

## Get info from design file
DESIGN_FILE=$PROJ/design_files/sampleID_ST_CC_ref_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
REF=${ARRAY[5]}

#### index bam files
samtools index ${ALN}/${SAMPLEID}_2_CCRef_unique_sorted.bam ${ALN}/${SAMPLEID}_2_CCRef_unique_sorted.bam.bai

## Call SNPs -- ploidy = 1 
freebayes -p 1 -f ${REF_DIR}/${REF}.fna -b ${ALN}/${SAMPLEID}_2_CCRef_unique_sorted.bam | vcffilter -f 'QUAL > 20' > ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}.vcf

## Rename sample columns to strain name
${PROJ}/scripts/vcf_rename_sample_01dbs.py -v ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}.vcf ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}_renamed.vcf ${SAMPLEID}
 
## Cull to SNPs only
vcftools --vcf ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}_renamed.vcf --remove-indels --recode --recode-INFO-all --out ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}_renamed_snps

## Cull to INDELS only
vcftools --vcf ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}_renamed.vcf --keep-only-indels --recode --recode-INFO-all --out ${OUTPUT}/${SAMPLEID}_2_CCRef_${REF}_renamed_indels
