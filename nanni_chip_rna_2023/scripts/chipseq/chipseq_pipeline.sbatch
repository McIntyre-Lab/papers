#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=chipseq
#SBATCH --output=/ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/chipseq_pipeline_%j.out
##SBATCH --time=72:00:00
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
##SBATCH --mem=16gb
#SBATCH --mem=4gb

date;hostname;pwd

### Set Global Variables

## Project, scripts, and temporary directories
export PROJ=/ufrc/mcintyre/share/etoh_srna
export SCRIPTS=$PROJ/scripts/chipseq
export ROZ=$PROJ/roz_chipseq_analysis
    mkdir -p ${ROZ}

## Initial Design File
export DESIGN_FILE=${PROJ}/design_files/design_etoh_srna_chipSeq_noHeader.csv

## Species Names (will be used in output file names)
export SPECIES1="mel"
export SPECIES2="sim"

## Corresponding FASTA and GTF files for each species
export FASTA1=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
export GTF1=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-r6.17.gtf
export FASTA2=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
export GTF2=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-r2.02.gtf

## Feature Annotation Directory
## (contains ${SPECIES}_*.unique.bed files for different genomic features)
#export FEATURES=$PROJ/ref_bed
#export FEATURES=/ufrc/mcintyre/share/etoh_srna/TEST_ANNOT_PIPELINE/ref_bed
#export FEATURES=/ufrc/mcintyre/share/etoh_srna/REDO_annotation_pipeline/ref_bed
export FEATURES=/ufrc/mcintyre/share/etoh_srna/dros_annotation_pipeline

## Corresponding fragment (exonic), intron, and intergenic BED files
export FRAGMENT1=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_annotations_100bp_reads/deml617_100bp_exon_fragments_coverage.bed
export FRAGMENT2=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_100bp_reads/dsim202_100bp_exon_fragments_coverage.bed
export INTRON1=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_annotations_100bp_reads/deml617_100bp_introns_from_fusions.bed
export INTRON2=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_100bp_reads/dsim202_100bp_introns_from_fusions.bed
export FUSION1=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_annotations_100bp_reads/deml617_100bp_fusions_coverage.bed
export FUSION2=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_100bp_reads/dsim202_100bp_fusions_coverage.bed
export INTERGEN1=${FEATURES}/${SPECIES1}_intergenic.bed
export INTERGEN2=${FEATURES}/${SPECIES2}_intergenic.bed

## Alignment Directory
export ALN=$PROJ/chipseq/aln_genome

#########################
###   Peak Analysis   ###
#########################
#source ${SCRIPTS}/chipseq_peak_analysis.sh


############################
###   Feature Analysis   ###
############################
source ${SCRIPTS}/chipseq_feature_analysis.sh


## Remove temporary directory
#rm -r ${ROZ}

