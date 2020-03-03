#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=refIndex
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/refIndex_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre



module purge
module load rsem/1.2.31 star/2.7.0b bwa/0.7.7 gcc/5.2.0 python/2.7.10 samtools/1.4 bowtie/0.12.9

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

OUT=$PROJ/pacbio_analysis/
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi

GTF=$PROJ/pacbio_analysis/merge_samples_b73/all_samples.chained.gff
FA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

### (1) create PB transcriptome using rsem-prepare-reference --> outputs fasta file
#echo "Prepare RSEM index... $(date)
#GTF used : ${GTF}
#REF used : $FA
#"
#rsem-prepare-reference -p 8 \
#    --gtf ${GTF} ${FA} $OUT/PB_rsem_index
#echo "prep ref done."

### (2) create PB transcriptome using rsem-prepare-reference for bowtie1
echo "Prepare RSEM BT1 index... $(date)
GTF used : ${GTF}
REF used : $FA
"
rsem-prepare-reference -p 8 \
    --gtf ${GTF} ${FA} $OUT/PB_rsem_index_BT1 --bowtie
echo "prep BT1 ref done."



#### (2) create bwa ref for PB transcriptome
##BWA-mem references  -p is the prefix for the output database, -a is the indexing algorithm ('bwtsw' is for ref>2G, 'is' for ref<2G).
#echo "Prepare BWA index from :  $OUT/PB_rsem_index.idx.fa
"
#bwa index -p $OUT/PB_rsem_index_BWA -a is  $OUT/PB_rsem_index.idx.fa


#### (3) samtools index for PB transcriptome
#echo "Preparing Samtools index "
#samtools faidx $OUT/PB_rsem_index.idx.fa


### (4) bowtie index for PB transcriptome
#echo "Preparing bowtie1 index"
#bowtie-build -f $OUT/PB_rsem_index.idx.fa $OUT/PB_rsem_index_BT1
