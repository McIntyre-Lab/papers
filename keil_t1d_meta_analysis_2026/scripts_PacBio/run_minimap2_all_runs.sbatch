#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=minimap2_allruns.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/minimap2_hq.%A_%a.out
#SBATCH--error=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/minimap2_hq.%A_%a.err
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --array=2-61
#SBATCH --account=concannon
#SBATCH --qos=concannon-b


module purge
module load minimap2/2.24
module load samtools/1.10


## align converted fq files to hg38 genome reference with minimap

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio

REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa
REF_GTF=/blue/mcintyre/share/transcript_distance/human_analysis/TranD_pairwise_ensembl_refseq_union/consolidated_transcriptome.gtf

## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_bc_2_run.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]} 
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}
SAMPLEID=${ARRAY[3]}
BC=${ARRAY[4]}

IND=$PROJ/isoseq3_refine_by_sample_isoseq/${RUN}/bam2fq

OUTD=${PROJ}/minimap2_all_runs_hq/${RUN}
    mkdir -p $OUTD
LOG=${OUTD}/log
    mkdir -p $LOG


#For mapping and prioritising annotated junctions
#Comment out because bed file only needs to be created once

#paftools.js gff2bed ${REF_GTF} > ${PROJ}/minimap2_all_runs_hq/minimap_hg38_union_ref.bed

#alignment using 
minimap2 \
    -t 12 \
    -ax splice:hq \
    --junc-bed ${PROJ}/minimap2_all_runs_hq/minimap_hg38_union_ref.bed \
    -u f \
    --secondary=no \
    ${REF_FA} \
    $IND/${SAMPLEID}_${RUN}.flnc.fastq \
    > $OUTD/${SAMPLEID}_${RUN}_unfiltered_minimap2.sam 2>$LOG/${SAMPLEID}_${RUN}_unfiltered_minimap2.log

## include header, include reads without flag 2052 (omit unmapped and supplementary)
samtools view \
    -h \
    -F 2052 \
    $OUTD/${SAMPLEID}_${RUN}_unfiltered_minimap2.sam \
    > $OUTD/${SAMPLEID}_${RUN}_minimap2.sam

rm $OUTD/${SAMPLEID}_${RUN}_unfiltered_minimap2.sam

