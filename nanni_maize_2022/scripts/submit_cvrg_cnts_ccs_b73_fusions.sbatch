#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=minimap2
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/cvr_cnts_%A_%a.out
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=12
#SBATCH --array=1-11

## fusion coverage counts of ccs reads mapped to B73 genome 


### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
ANNO=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73_EA_150bp
ALN=$PROJ/check_isoseq3_lost_genes/minimap2

BAMFILES=${ALN}/bamfiles
    mkdir -p $BAMFILES
SAMFILES=${ALN}/samfiles
    mkdir -p $SAMFILES
PILES=${ALN}/mpileups
    mkdir -p $PILES
CVR=$PROJ/cvr_cnts_ccs_b73_fusions
    mkdir -p $CVR
LOGC=$CVR/log_files
    mkdir -p $LOGC

## Get reference genome FASTA
FASTA=/blue/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa
BED=${ANNO}/FSM_consolidation_maize_B73_fusions_coverage.bed

## Get info from design file
DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}

date

## try coverage counts

    ROZ=${PROJ}/roz_${SAMPLE_NAME}
        mkdir -p ${ROZ}

    module load bwa/0.7.7 gcc/5.2.0 python/2.7.10 samtools

    ## sam to bam
    samtools view -ut ${FASTA}.fai -o ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.bam ${ALN}/${SAMPLE_NAME}_minimap2_b73.sam

    ## create sorted bam and index
    samtools sort -T ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.tmp.sorted \
        -o ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.bam \
        -O bam ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.bam
    samtools index ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.bam

    ## create mpileup from bam files
    samtools index ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.bam
    samtools mpileup -d 1000000 -f ${FASTA} ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.bam \
        > ${PILES}/${SAMPLE_NAME}_minimap2_b73.sorted.mpileup

    ## convert sorted bam to sorted sam
    samtools view -h -o ${SAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.sam ${BAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.bam

    ## get coverage counts for alignments using rpkm_calculate
    python ${PROJ}/scripts/rpkm_calculate.py \
        -m ${PILES}/${SAMPLE_NAME}_minimap2_b73.sorted.mpileup \
        -s ${SAMFILES}/${SAMPLE_NAME}_minimap2_b73.sorted.sam \
        -b ${BED} \
        -g ${LOGC}/cvrg_cnts_${SAMPLE_NAME}.logfile \
        -o ${ROZ}/cvrg_cnts_${SAMPLE_NAME}.csv \
        -c

    ## add column with sampleID to cov counts
    awk -F "," -v id="${SAMPLE_NAME}" -v OFS="," '{print id, $0}' ${ROZ}/cvrg_cnts_${SAMPLE_NAME}.csv > ${CVR}/cvrg_cnts_${SAMPLE_NAME}.csv


done
