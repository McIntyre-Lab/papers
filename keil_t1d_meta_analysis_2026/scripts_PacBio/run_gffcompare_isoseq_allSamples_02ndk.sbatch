#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=gff_comp2.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/gff_comp_all_genes.%A_%a.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/gff_comp_all_genes.%A_%a.err 
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=${PROJ}/isoseq3cluster_allSamples/sam2gtf
OUTD=${PROJ}/gff_compare_isoseq_allSamples_02ndk
    mkdir -p ${OUTD}

module load gffcompare/0.12.2

REF=/blue/mcintyre/share/transcript_distance/human_analysis/TranD_pairwise_ensembl_refseq_union/consolidated_transcriptome.gtf

echo "start gff compare"
date

gffcompare -r ${REF} -o ${OUTD}/isoseq_allSamples  ${IND}/allSamples_realign.gtf

echo "gffcompare done"
date

