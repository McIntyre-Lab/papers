#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=df
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/df_%A_%a.out
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1

## create design file with 1 line per unpolished bam file preceded by ID, geno, and trt

PROJ=/ufrc/mcintyre/share/maize_ainsworth

#for FILE in $PROJ/isoseq3_analysis/cluster_by_individual/*; do find $FILE \( -name "*?.bam" -o -name "*.??.bam" \); done | awk 'BEGIN{FS="/"}{split($NF,path,"_"); print path[1]","path[2]","substr(path[3],1,3)","$NF}' | sort > $PROJ/design_files/design_maize_PacBio_cluster_test_unpolished_path_all.csv

## remove rows containing 'flnc'
##sed -i '/flnc/d' $PROJ/design_files/design_maize_PacBio_cluster_test_unpolished_path_all.csv

