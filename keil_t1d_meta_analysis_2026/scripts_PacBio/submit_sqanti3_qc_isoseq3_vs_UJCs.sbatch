#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti_UJCs
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sqqanti_UJC.%j.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=1-25

## run sqanti3 QC

module load sqanti3/5.1.1

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
UJC_DIR=${PROJ}/id_UJC

## design does not have header!!
DESIGN_FILE=$PROJ/design_files/design_chromosomes.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

CHROM=${ARRAY[0]}


REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa
GTF=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk/isoseq_allSamples_corrected_geneID.gtf

OUTD=${PROJ}/sqanti3_qc_isoseq_vs_UJCs
    mkdir -p $OUTD
ROZ=${OUTD}/roz
    mkdir -p ${ROZ}


##Subset isoseq3 gtf by chromosome

awk -v chrom=${CHROM} '{if($1==chrom) {print $0}}' ${GTF} >  ${ROZ}/chrom${CHROM}_isoseq_allSamples_corrected_geneID.gtf

##Subset UJC gtf by chromosome

awk -v chrom=${CHROM} '{if($1==chrom) {print $0}}' ${UJC_DIR}/allSamples_UJC.gtf >  ${ROZ}/chrom${CHROM}_allSamples_UJC.gtf

##RUN SQANTI AGAINST UJCs GTF

##  --gtf option indicates GTF isoform input
##  --skipORF option to skip ORF prediction (runs much faster)
##  --force_id_ignore option to allow transcript ID names that are not like PB.X.Y
##  Order of files is (1) isoforms (FA/FQ/GTF), (2) annotation (ref GTF/GFF), (3) genome (ref FA)

#   --dir is path to output directory
#   -o prefix name to use

## order of files input into sqanti3_qc.py is critical!!
# 1. sam2gtf (from read align to genome)
# 2. FSM consolidated ref gtf file
# 3. ref genome fasta file


sqanti3_qc.py \
    --skipORF \
    --force_id_ignore \
    --dir ${OUTD} \
    -o ${CHROM} \
    ${ROZ}/chrom${CHROM}_isoseq_allSamples_corrected_geneID.gtf \
    ${ROZ}/chrom${CHROM}_allSamples_UJC.gtf \
    ${REF_FA}


rm ${ROZ}/chrom${CHROM}_isoseq_allSamples_corrected_geneID.gtf
rm ${ROZ}/chrom${CHROM}_allSamples_UJC.gtf
