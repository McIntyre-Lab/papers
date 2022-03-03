#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=prep_rsem
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/prep_rsem_index_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcintyre

module purge
module load rsem/1.2.31 star/2.7.5c bwa/0.7.7 gcc/5.2.0 bowtie/0.12.9

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth

OUT=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/references
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi

GTF=$PROJ/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf
FA=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

### (1) create reference genome index using rsem-prepare-reference
echo "Prepare RSEM index... $(date)
GTF used : ${GTF}
REF used : ${FA}
"
rsem-prepare-reference -p 8 \
    --gtf ${GTF} ${FA} $OUT/fsm_ism_nic_nnc_consol_junc_rsem_index
echo "prep ref done."

### (2) create reference genome index using rsem-prepare-reference for star
echo "Prepare RSEM STAR index... $(date)
GTF used : ${GTF}
REF used : $FA
"
rsem-prepare-reference -p 8 \
    --gtf ${GTF} ${FA} $OUT/fsm_ism_nic_nnc_consol_junc_rsem_index_STAR --star
echo "prep STAR ref done."
