#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sq_merge_mo17
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/sq_merge_mo17_%A.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1

### SQANTI of merged maize PacBio data - mapped to Mo17_Yan

module load pacbio/5.1.0 isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

REF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.fasta

REF_GTF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.mw2.gtf

## Set input and output directories
OUT=$PROJ/sqanti_merged_mo17_yan
    mkdir $OUT
IN=$PROJ/merge_samples_mo17_yan

## GTF (is GTF but has .gff ending) from output of collapse_isoforms_by_sam.py
PACBIO_GTF=$IN/all_samples.chained.gff

date
echo "***SQANTI***
"
python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py -g ${PACBIO_GTF} ${REF_GTF} \
    $REF -o maize -d $OUT 2>$OUT/maize.sqanti.log
