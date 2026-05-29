#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=Isoseq3Cluster_allSamples.%A_%a
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/isoseq3Cluster_allSamples.%j.out
#SBATCH --account=concannon
#SBATCH --qos=concannon
#SBATCH --cpus-per-task=24
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --partition=bigmem
#SBATCH --mem=400gb

module load isoseq3/3.4.0
module load samtools/1.10
module load minimap2/2.24
module load bedtools/2.29.2

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=$PROJ/baseline_selection_keep_genic
OUTD=$PROJ/isoseq3cluster_allSamples
    mkdir -p ${OUTD}
ROZ=${OUTD}/roz
    mkdir -p ${ROZ} 

LOG=${OUTD}/log
    mkdir -p ${LOG}

ALIGN_OUTD=$OUTD/minimap2align
ALIGN_LOG=${ALIGN_OUTD}/log
    mkdir -p ${ALIGN_LOG}

sam2gtf_OUTD=$OUTD/sam2gtf
    mkdir -p ${sam2gtf_OUTD}


REFDIR=/blue/concannon/share/jnewman/references/hg38_ensembl
REF=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa


echo "
*** Combine all samples ***
"
samtools merge ${ROZ}/temp.bam \
    ${IND}/*RNA_baseline_select_keep_genic.bam

samtools merge ${ROZ}/allSamples.bam \
    ${ROZ}/temp.bam ${IND}/*case_baseline_select_keep_genic.bam

    echo "
*** IsoSeq3 transcript estimation ***
"
##make fofn file
find ${IND}/*RNA_baseline_select_keep_genic.bam -type f > ${ROZ}/allSamples.fofn
find ${IND}/*case_baseline_select_keep_genic.bam -type f >> ${ROZ}/allSamples.fofn

### run isoseq3 cluster
isoseq3 cluster \
    --log-file ${LOG}/allSamples.log \
    --verbose \
    --use-qvs \
    ${ROZ}/allSamples.bam \
    ${OUTD}/allSamples_isoseq3_clustered.bam

echo "
*** Map isoseq transcripts using minimap2 ***
"

minimap2 \
    -t 12 \
    -a \
    -x splice \
    -u f \
    --secondary=no \
    -C 5 \
    ${REF} \
    ${OUTD}/allSamples_isoseq3_clustered.lq.fasta.gz \
    > ${ALIGN_OUTD}/allSamples_realigned.lq.sam 2>${ALIGN_LOG}/allSamples_realigned.lq.log

minimap2 \
    -t 12 \
    -a \
    -x splice \
    -u f \
    --secondary=no \
    -C 5 \
    ${REF} \
    ${OUTD}/allSamples_isoseq3_clustered.hq.fasta.gz \
    > ${ALIGN_OUTD}/allSamples_realigned.hq.sam 2>${ALIGN_LOG}/allSamples_realigned.hq.log

#### sam to gtf
samtools merge \
    ${sam2gtf_OUTD}/allSamples_realigned.sam \
    ${ALIGN_OUTD}/allSamples_realigned.hq.sam \
    ${ALIGN_OUTD}/allSamples_realigned.lq.sam
samtools view \
    -b \
    -h \
    -F 2052 \
    ${sam2gtf_OUTD}/allSamples_realigned.sam \
    > ${sam2gtf_OUTD}/allSamples_realigned.bam

bedtools bamtobed \
    -split \
    -i ${sam2gtf_OUTD}/allSamples_realigned.bam \
    | awk -F "\t" -v source="allSamples" \
    '{print $1"\t"source"\t""exon""\t"$2+1"\t"$3"\t"".""\t"$6"\t"".""\t""gene_id " "\"" $4 "\"" "; " "transcript_id " "\"" $4 "\"" ";"}' \
    > ${sam2gtf_OUTD}/allSamples_realign.gtf

rm ${sam2gtf_OUTD}/allSamples_realigned.sam

