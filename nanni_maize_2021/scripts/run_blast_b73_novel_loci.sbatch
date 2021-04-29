#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/blast_b73_only_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## BLAST B73 exclusive transcripts from novel loci with at least
##     5 reads in both B73 treatments

module purge
module load ncbi_blast/2.10.1

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INDIR=$PROJ/htseq_quantify_genes_loci

blastn \
    -db nt \
    -query ${INDIR}/novel_loci_htseq_B73_only_novel_loci_ge5reads.fasta \
    -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length bitscore evalue mismatch gaps" \
    > ${INDIR}/blast_novel_lovi_B73_only_ge5reads.tsv

awk -F "\t" '{if(name!=$1){print $0; name=$1}}' \
    ${INDIR}/blast_novel_lovi_B73_only_ge5reads.tsv \
    ${INDIR}/blast_novel_lovi_B73_only_ge5reads_top_hit.tsv
