#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=target
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/staph_relapse/scripts/SLURM_LOGS/target_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb    		# per processor memory


## create list of snps for each reference


## Set directories
PROJ=/ufrc/mcintyre/share/staph_relapse
REF_DIR=${PROJ}/CC_references
VCF_IN=${PROJ}/vcf_by_ref
BAM=$PROJ/aln_bwa_PE


module load gcc samtools/1.4 bcftools bwa/0.7.7 python/2.7.6

for REF in MSSA476 ED98 CA_347 MSSA476 Newman ST20130941 TCH60
do

    ## create targets file 
    bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${VCF_IN}/${REF}_snps.vcf | bgzip -c > ${VCF_IN}/${REF}_snps.tsv.gz && tabix -s1 -b2 -e2 ${VCF_IN}/${REF}_snps.tsv.gz

done
