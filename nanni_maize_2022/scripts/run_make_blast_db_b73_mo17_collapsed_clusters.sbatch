#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast_db
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/blast_db_clusters_%A.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

### make blast db for b73 and mo17 collapsed clusters

module load ncbi_blast seqtk

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INPUT=$PROJ/compare_b73_2_mo17
OUTPUT=${INPUT}/blast_db_collapsed_clusters
    mkdir -p $OUTPUT

## Set up temp directory
ROZ=${OUTPUT}/roz_blast_db
    mkdir -p ${ROZ}

for REF in b73; do
#for REF in b73 mo17_yan mo17_cau; do
    for COMPARE in b73 mo17; do
        ## Convert collapsed cluster FQ to FA
        seqtk seq -A \
            ${INPUT}/tofu2_${REF}/*${COMPARE}*/*.collapsed.rep.fq \
            > ${ROZ}/${COMPARE}_collapsed.fa
        makeblastdb -in ${ROZ}/${COMPARE}_collapsed.fa \
            -dbtype 'nucl' \
            -out ${OUTPUT}/${COMPARE}_collapsed.fa.BLAST
    done
done

rm -r ${ROZ}
