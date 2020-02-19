#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=minimap2
#SBATCH	--output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/map_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11

### Minimap2 mapping of maize PacBio data 
###  Aln to b73 ref genome fasta sequence (DNA seq, top level, soft masked)

### all samples

module load minimap/2.12

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTD=$PROJ/mapping_minimap2_b73
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IND=$PROJ/isoseq3_analysis/polish_by_individual

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

## Set input and output directories to the correct species
OUT=$OUTD/$GENO/$TRT
    if [ ! -e $OUT ]; then mkdir -p $OUT; fi
IN=$IND/$GENO/$TRT

date
echo "***Combine Polish Fasta.gz***
"
## Combine Polish hq.fasta.gz output files
## hq files are the high quality clusted consensus reads with predicted accuracy >= 0.99

if [ ! -e $OUT/${SAMPLE_NAME}.polished.all.hq.fasta.gz ]; then
    zcat $IN/${SAMPLE_NAME}.*.?.hq.fasta.gz $IN/${SAMPLE_NAME}.*.??.hq.fasta.gz > $OUT/${SAMPLE_NAME}.polished.all.hq.fasta
    echo "    Creating combined polish fasta file"
else
    echo "    Combined polish fasta file already exists"
fi

date
echo "***Combine Polish Fastq.gz***
"
## Combine Polish hq.fastq.gz output files
## Not used in mapping but used in downstream analysis
## hq files are the high quality clusted consensus reads with predicted accuracy >= 0.99

if [ ! -e $OUT/${SAMPLE_NAME}.polished.all.hq.fastq.gz ]; then
    zcat $IN/${SAMPLE_NAME}.*.?.hq.fastq.gz $IN/${SAMPLE_NAME}.*.??.hq.fastq.gz > $OUT/${SAMPLE_NAME}.polished.all.hq.fastq
    echo "    Creating combined polish fastq file"
else
    echo "    Combined polish fastq file already exists"
fi

date
echo "***Minimap2 Map***
"
## Minimap2 Map
## -a output in SAM format
## -x splice uses long read spliced alignment settings
##     meaning long deletions are considered introns, no long insertions,
##     and gap costs are changed accordingly
## -u f option requires mapping of the forward strand
## --secondary=no does not include secondary alignments
## -C 5 modify cost of non-canonical splice junctions (recommended for Isoseq processed PacBio)

minimap2 -t 12 -a -x splice -u f --secondary=no -C 5 $REF $OUT/${SAMPLE_NAME}.polished.all.hq.fasta \
    > $OUT/${SAMPLE_NAME}.polished.all.hq.mapped.sam 2>$OUT/${SAMPLE_NAME}.minimap2.log

date
