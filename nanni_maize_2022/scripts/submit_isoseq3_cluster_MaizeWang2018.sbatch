#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=iso_clus
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/iso_lib_proc_MaizeWang2018_cluster_%A_%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=12
#SBATCH --array=2-9

### isoseq3 analysis of MaizeWang2018 samples (pool 1) - library processing
### 9 samples total

module purge
module load samtools/1.10 isoseq3/3.0.0

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/isoseq3_library_processing
LIMAOUT=${OUTD}/lima
CLUSTEROUT=${OUTD}/cluster
    mkdir -p ${CLUSTEROUT}

## Get info from design file

DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Make temporary directory for concatenated bam files
ROZ=${CLUSTEROUT}/roz_${SAMPLEID}
    mkdir -p ${ROZ}

date

echo "
***Cluster By Individual***
"
## Concatenate all demultiplexed bam from the same sample
samtools merge -f ${ROZ}/${SAMPLEID}_ccs_demux.primer_5p--${BC}_3p.bam \
    ${LIMAOUT}/*_ccs_demux.primer_5p--${BC}_3p.bam

## Cluster
## input is the consensus read bam *.ccs.bam
## removes polyA and artificial concatemers to create FLNC (Full-Length,Non-Concatemer) reads, then puts them in transcript sets
## --split-bam 24 splits the output into 24 files for downstream parallelization
## outputs unpolished.#.bam and unpolished.transcriptset.#.xml where # is the file number after splitting
## use --require-polya --> trims polyA seq

isoseq3 cluster \
    --require-polya \
    --num-threads 12 \
    -v \
    --split-bam 24 \
    ${ROZ}/${SAMPLEID}_ccs_demux.primer_5p--${BC}_3p.bam \
    ${CLUSTEROUT}/${SAMPLEID}.unpolished.bam

date

rm -r ${ROZ}
