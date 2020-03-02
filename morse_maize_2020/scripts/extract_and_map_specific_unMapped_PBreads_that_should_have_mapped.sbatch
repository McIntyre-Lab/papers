#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=extract
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/extract_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --cpus-per-task=12
#SBATCH --account=mcintyre


## extract transcript/9850 from 113_c123 amb with mo17 yan ref (this looks like and excellent blast hit to cds)
## map to all 3 references 

module load samtools minimap2.12

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
INPUT=$PROJ/uncollapsed_mo17_yan_IsoSeq_reads
OUTPUT=$PROJ/minimap2_should_haves 
    mkdir -p $OUTPUT

## extract sequence
#samtools faidx $INPUT/113_c123_amb_unCollapsed_unmapped_mo17_yan_ref.fa
#samtools faidx $INPUT/113_c123_amb_unCollapsed_unmapped_mo17_yan_ref.fa transcript/9850 \
#    > $INPUT/transcript_9850_113_C123_unmapped_mo17_yan_ref.fa


## Minimap2 Map
## -a output in SAM format
## -x splice uses long read spliced alignment settings
##     meaning long deletions are considered introns, no long insertions,
##     and gap costs are changed accordingly
## -x asm5 long assembly to reference mapping.  typically alignment will not extend to regions with 5% or higher seq divergence
## -x asm10 long assembly to reference mapping. up to 10% seq divergence
## -x asm20 long assembly to reference mapping. up to 20% seq divergence
## -u f option requires mapping of the forward strand
## --secondary=no does not include secondary alignments
## -C 5 modify cost of non-canonical splice junctions (recommended for Isoseq processed PacBio)
## --paf-no-hit output unmapped queries in PAF

B73=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
MO17_YAN=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.fasta
MO17_CAU=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa

echo "query     query_length    query_start     query_end	strand  target  target_length target_start	target_end	num_match	num_mapped	mapping_qual    tp	cm     $

## (1) map to B73 with splice alignment settings
minimap2 \
    -t 12 \
    -x splice \
    -u f \
    --secondary=yes \
    -C 5 \
    $B73 \
    $INPUT/transcript_9850_113_C123_unmapped_mo17_yan_ref.fa \
    > $ROZ/B73_splice_mapped.paf 2>$OUTPUT/B73_splice_mapped.log

cat $ROZ/header.txt $ROZ/B73_splice_mapped.paf > $OUTPUT/B73_splice_mapped.paf

## (2) map to Mo17 Yan with splice alignment settings
minimap2 \
    -t 12 \
    -x splice \
    -u f \
    --secondary=yes \
    -C 5 \
    $MO17_YAN \
    $INPUT/transcript_9850_113_C123_unmapped_mo17_yan_ref.fa \
    > $ROZ/Mo17_Yan_splice_mapped.paf 2>$OUTPUT/Mo17_Yan_splice_mapped.log

cat $ROZ/header.txt $ROZ/Mo17_Yan_splice_mapped.paf > $OUTPUT/Mo17_Yan_splice_mapped.paf

## (3) map to Mo17 Cau with splice alignment settings
minimap2 \
    -t 12 \
    -x splice \
    -u f \
    --secondary=yes \
    -C 5 \
    $MO17_CAU \
    $INPUT/transcript_9850_113_C123_unmapped_mo17_yan_ref.fa \
    > $ROZ/Mo17_Cau_splice_mapped.paf 2>$OUTPUT/Mo17_Cau_splice_mapped.log

cat $ROZ/header.txt $ROZ/Mo17_Cau_splice_mapped.paf > $OUTPUT/Mo17_Cau_splice_mapped.paf



