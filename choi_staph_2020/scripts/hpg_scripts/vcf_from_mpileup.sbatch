#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=vcf
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/staph_relapse/scripts/SLURM_LOGS/vcf_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb    		# per processor memory
#SBATCH --array=1-78

### convert sorted bam to vcf for snps in targets file

## Set directories
PROJ=/ufrc/mcintyre/share/staph_relapse
REF_DIR=${PROJ}/CC_references
VCF_IN=${PROJ}/vcf_by_ref
BAM=$PROJ/aln_bwa_PE

## Get info from design file
DESIGN_FILE=${PROJ}/design_files/sampleID_ST_CC_ref_noS15pair_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
SAMPLENUM=${ARRAY[1]}
REF=${ARRAY[5]}
TAXONOMYID=${ARRAY[6]}


module load gcc samtools/1.4 bcftools bwa/0.7.7 python/2.7.6

VCF=${PROJ}/vcf_w_targets_file
    mkdir -p ${VCF}

## Convert sorted BAM to VCF for snps in target file
#samtools mpileup -d 10000000000 -uf ${REF_DIR}/${REF}.fna ${BAM}/${SAMPLEID}_2_CCRef_unique_sorted.bam | bcftools call -m --ploidy 1 -T ${VCF_IN}/${REF}_snps.tsv.gz - > ${VCF}/${SAMPLEID}_2_CCRef_snps.vcf

## rename sample columns in vcf to sampleID
#${PROJ}/scripts/vcf_rename_sample_01dbs.py -v ${VCF}/${SAMPLEID}_2_CCRef_snps.vcf ${VCF}/${SAMPLEID}_2_CCRef_${REF}_renamed_snps.vcf ${SAMPLEID}

sed '/^##/d' ${VCF}/${SAMPLEID}_2_CCRef_${REF}_renamed_snps.vcf ${SAMPLEID} > ${VCF}/${SAMPLEID}_2_CCRef_${REF}_renamed_snps_noHeader.vcf ${SAMPLEID}
