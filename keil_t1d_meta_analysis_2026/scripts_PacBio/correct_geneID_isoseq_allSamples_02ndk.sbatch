#!/bin/bash
#SBATCH--mail-user=nkeil@ufl.edu
#SBATCH--mail-type=FAIL,END
#SBATCH--job-name=correct_gene.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/correct_allgenes.%A_%a.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/correct_allgenes.%A_%a.err 
#SBATCH--cpus-per-task=1
#SBATCH--time=24:00:00
#SBATCH--ntasks=1
#SBATCH--mem=256gb
#SBATCH--account=concannon
#SBATCH--qos=concannon-b



### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
GFFCOMP=${PROJ}/gff_compare_isoseq_allSamples_02ndk
IND=${PROJ}/isoseq3cluster_allSamples/sam2gtf

OUTD=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk
    mkdir -p $OUTD

module purge
module load python/3.6


#Correct all samples GTF
echo "starting correction"
date

## Correct gene_id values using associated gene_id or xloc id from gffcompare
python $PROJ/scripts/correct_gffcompare_GTF_gene_id_06avn.py \
    -a ${GFFCOMP}/isoseq_allSamples.annotated.gtf \
    -g ${IND}/allSamples_realign.gtf \
    --key ${OUTD}/isoseq_gene_id_2_gffcomp_gene_id.csv \
    -o ${OUTD}/isoseq_allSamples_corrected_geneID.gtf

echo "correction done"
date

