#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=merge
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/merge_%A_%a.out
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1

### Use chain_samples.py to combine all samples for comparison
### Must first summarize junctions and "scrub" them together


module load python tofu2/20180629

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth

REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

## Set output directory
OUT=$PROJ/merge_samples_b73
    mkdir -p $OUT

date

## Copy all required tofu2 output files to another directory - names need to be identical for combining
for P1 in $PROJ/tofu2_b73/*
do
    echo $P1
    ID=$(basename $P1)
    echo $ID
    OUT1=$PROJ/merge_samples_b73/config/${ID}
        mkdir -p $OUT1
    cp "$PROJ/tofu2_b73/${ID}/${ID}.collapsed.gff" $OUT1/original.gff
    cp "$PROJ/tofu2_b73/${ID}/${ID}.collapsed.group.txt" $OUT1/original.group.txt
    cp "$PROJ/tofu2_b73/${ID}/${ID}.collapsed.rep.fq" $OUT1/original.rep.fq
    cp "$PROJ/tofu2_b73/${ID}/${ID}.collapsed.abundance.txt" $OUT1/original.count.txt
done
 
## Make summarize.config file
for P1 in $PROJ/merge_samples_b73/config/*
do 
    echo "SAMPLE=$(basename $P1);$P1/"
done > $OUT/maize_summarize.config
echo "GFF_FILENAME=original.gff" >> $OUT/maize_summarize.config
echo "GENOME_FILENAME=${REF}" >> $OUT/maize_summarize.config


## Summarize Junctions of all tissues
summarize_sample_GFF_junctions.py $OUT/maize_summarize.config \
    $OUT/maize

## Modify summarized file so that each junction is all on one line
## Default puts label on a separate line for some reason
awk '{if($1=="chr"){print;} else if($2){printf "%s+",$0;} else{print $0}}' \
    $OUT/maize.junction_detail.txt > \
    $OUT/maize.junction_detail.mod.txt

## Make scrub.config file
for path in $PROJ/merge_samples_b73/config/* ; do
    echo "SAMPLE=$(basename $path);$path/"
done > $OUT/maize_scrub.config
echo "GFF_FILENAME=original.gff" >> $OUT/maize_scrub.config
echo "COUNT_FILENAME=original.count.txt" >> $OUT/maize_scrub.config
echo "GROUP_FILENAME=original.group.txt" >> $OUT/maize_scrub.config
echo "FASTQ_FILENAME=original.rep.fq" >> $OUT/maize_scrub.config

## Make temp directory and mv into it
## Outputs are the same and would overwrite each other
mkdir  temp_maize_b73
cd temp_maize_b73

## Scrub junctions of each sample
scrub_sample_GFF_junctions.py $OUT/maize_scrub.config \
    $OUT/maize.junction_detail.mod.txt scrubbed

## Output structure requires there be no path for output prefix
## Must move the one output file to $OUT
mv scrubbed.scrubbed.junction.bed $OUT/maize_scrubbed.junction.bed

date
echo "
*** Chaining***
"

## Make sample.config file
for path in $PROJ/merge_samples_b73/config/* ; do
    echo "SAMPLE=$(basename $path);$path/"
done > $OUT/maize_sample.config
echo "GFF_FILENAME=scrubbed.gff" >> $OUT/maize_sample.config
echo "COUNT_FILENAME=scrubbed.abundance.txt" >> $OUT/maize_sample.config
echo "GROUP_FILENAME=scrubbed.group.txt" >> $OUT/maize_sample.config
echo "FASTQ_FILENAME=scrubbed.rep.fq" >> $OUT/maize_sample.config


## Chain tissue samples together for each organism
## Can use the same config file from scrubbing
## creates its own output files all_samples.chained
chain_samples.py $OUT/maize_sample.config norm_nfl

## Move files out of temp directory and remove it
mv * $OUT/
cd ..
rm -r temp_maize_b73

date

