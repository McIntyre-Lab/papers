#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=aln_2_iso
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/staph_relapse/scripts/SLURM_LOGS/aln_2_%A-%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb    		# per processor memory
#SBATCH --array=1-79


#### align (bwa mem PE) isolate reads to references based on CC
#### sam to sorted bam

## Set directories
PROJ=/ufrc/mcintyre/share/staph_relapse
REF_DIR=${PROJ}/CC_references

## Get info from design file
DESIGN_FILE=${PROJ}/design_files/sampleID_ST_CC_ref_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
SAMPLENUM=${ARRAY[1]}
REF=${ARRAY[5]}
TAXONOMYID=${ARRAY[6]}

## Create output directory
OUTPUT=${PROJ}/aln_bwa_PE
    mkdir -p $OUTPUT
ROZ=${PROJ}/roz_${SAMPLEID}
    mkdir -p $ROZ

module load gcc samtools/1.4 bwa/0.7.7 python/2.7.6

READS=${PROJ}/fastq_split_dups_PE

## BWA-mem alignment
#bwa mem -t 8 ${PROJ}/CC_references/${REF}.fna_bwa ${READS}/*-${SAMPLENUM}_*R1_001_distinct.fq ${READS}/*-${SAMPLENUM}_*_R2_001_distinct.fq > ${OUTPUT}/${SAMPLEID}_2_CCRef.sam

## Split SAM file
#python ${PROJ}/scripts/BWASplitSAM_05dbs.py \
#    -s ${OUTPUT}/${SAMPLEID}_2_CCRef.sam \
#    -fq1 ${READS}/*-${SAMPLENUM}_*_L002_R1_001_distinct.fq \
#    -fq2 ${READS}/*-${SAMPLENUM}_*_L002_R1_001_distinct.fq \
#    --outdir ${OUTPUT}

## Concatenate mapped and opposite
#cat ${OUTPUT}/${SAMPLEID}_2_CCRef_mapped.sam ${OUTPUT}/${SAMPLEID}_2_CCRef_opposite.sam > ${OUTPUT}/${SAMPLEID}_2_CCRef_unique.sam

## Convert SAM to BAM
#samtools view -bt ${REF_DIR}/${REF}.fna.fai -o ${OUTPUT}/${SAMPLEID}_2_CCRef_unique.bam ${OUTPUT}/${SAMPLEID}_2_CCRef_unique.sam

## Sort BAM and index
samtools sort -@ 8 -m 120M -o ${OUTPUT}/${SAMPLEID}_2_CCRef_unique_sorted.bam  -O bam ${OUTPUT}/${SAMPLEID}_2_CCRef_unique.bam

## Remove all sam files for space reasons
rm ${OUTPUT}/${SAMPLEID}_2_CCRef*.sam

rm -r $ROZ
