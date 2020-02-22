#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/sqanti_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1

### SQANTI of merged maize PacBio data - mapped to B73

module load pacbio/5.1.0 isoseq3/3.0.0 samtools/1.3.1 gcc/7.3.0 R/3.2.0

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

REF_GTF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf

## Set input and output directories
OUT=$PROJ/sqanti_merged_b73
    mkdir $OUT
IN=$PROJ/merge_samples_b73

## GTF (is GTF but has .gff ending) from output of collapse_isoforms_by_sam.py
PACBIO_GTF=$IN/all_samples.chained.gff

date
echo "***SQANTI***
"
python $PROJ/scripts/pacbio/ConesaLab-sqanti/sqanti_qc.py -g ${PACBIO_GTF} ${REF_GTF} \
    $REF -o maize -d $OUT 2>$OUT/maize.sqanti.log
