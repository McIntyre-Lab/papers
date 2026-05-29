#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti_reads
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sqqanti_reads02.%j.out
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=400gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b

## run sqanti3 QC

module load sqanti3/5.1.1

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
READS_DIR=${PROJ}/consol_allSamples



REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa
GTF=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk/isoseq_allSamples_corrected_geneID.gtf

OUTD=${PROJ}/sqanti3_qc_isoseq_vs_reads_02ndk
    mkdir -p $OUTD

##Cat consolidated reads gtf together 

cat ${READS_DIR}/*/*.gtf > ${READS_DIR}/allSamples_allChr_consolidated_transcriptome.gtf

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
    -o allSamples_isoseq \
    ${GTF} \
    ${READS_DIR}/allSamples_allChr_consolidated_transcriptome.gtf \
    ${REF_FA}


