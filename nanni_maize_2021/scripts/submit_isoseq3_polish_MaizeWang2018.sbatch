#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_clus
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/iso_lib_proc_MaizeWang2018_polish_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=12
#SBATCH --array=2-216

### 216 unpolished cluster files
### isoseq3 analysis of MaizeWang2018 samples (pool 1) - library processing
### 9 samples total

module purge
module load isoseq3/3.0.0 pacbio/5.1.0

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/isoseq3_library_processing
BAMOUT=${OUTD}/subreads_bam
CLUSTEROUT=${OUTD}/cluster
POLISHOUT=${OUTD}/polish
    mkdir -p ${POLISHOUT}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_pool1_PacBio_cluster_unpolished_path.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SPECIES=${ARRAY[0]}
TISSUE=${ARRAY[1]}
UNPOLISHED=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Get unpolished file number
NUM=$(echo "${UNPOLISHED}" | awk -F "." '{print $(NF-1)}')

date

echo "
***Polish***
"

## Polish
## creates a consensus sequence for each clustered transcripts
## input is a unpolished.bam or a split (from isoseq3 cluster --split-bam) unpolished.#.bam
## outputs polished.transcriptset.#.bam

isoseq3 polish \
    ${CLUSTEROUT}/${UNPOLISHED} \
    ${BAMOUT}/MaizeWang2018_pool1_subreadset.xml \
    ${POLISHOUT}/${SAMPLEID}.polished.transcriptset.${NUM}.bam

date
