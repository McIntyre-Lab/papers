#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=prep_polish
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/prep_4_polish_MaizeWang2018_%A.out
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=1

## Create a design file with one line per unpolished bam file
##     with columns: species, tissue, bamfile
## Generate subread set XML file (using dataset create)
##     to be used as input for polishing step (generally is a
##     link to all subread.bam files)

module purge
module load isoseq3/3.0.0 pacbio/5.1.0

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/isoseq3_library_processing
BAMOUT=${OUTD}/subreads_bam
CLUSTEROUT=${OUTD}/cluster

## Make design file
for FILE in ${CLUSTEROUT}/*; do
    find ${FILE} \( -name "*.?.bam" -o -name "*.??.bam" \)
done | awk 'BEGIN{FS="/"}{split($NF,path,"_"); split(path[2],name,"."); \
    print path[1]","name[1]","$NF}' | sort \
    > $PROJ/design_files/df_MaizeWang2018_pool1_PacBio_cluster_unpolished_path.csv

## Make subread set
dataset create \
    --type SubreadSet \
    ${BAMOUT}/MaizeWang2018_pool1_subreadset.xml \
    ${BAMOUT}/*subreads.bam
