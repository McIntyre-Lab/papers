#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=refRsemStarI
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/refRsemStarI_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre



module purge
module load rsem/1.2.31 star/2.7.0b bwa/0.7.7 gcc/5.2.0 python/2.7.10 samtools/1.4 bowtie/0.12.9

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

OUT=$PROJ/pacbio_analysis/references
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi

GTF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
FA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

### (1) create reference genome index using rsem-prepare-reference
echo "Prepare RSEM index... $(date)
GTF used : ${GTF}
REF used : $FA
"
rsem-prepare-reference -p 8 \
    --gtf ${GTF} ${FA} $OUT/genome_rsem_index
echo "prep ref done."

### (2) create reference genome index using rsem-prepare-reference for star
echo "Prepare RSEM STAR index... $(date)
GTF used : ${GTF}
REF used : $FA
"
rsem-prepare-reference -p 8 \
    --gtf ${GTF} ${FA} $OUT/genome_rsem_index_STAR --star
echo "prep STAR ref done."
