#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --job-name=run_miso_bwa
#SBATCH --partition=hpg1-compute
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.aln_bwa.%j.%A.%a.out
#SBATCH -t 16:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30gb
#SBATCH --array=1-4
#

## Modules

module load miso/0.5.3
module load samtools
module load bedtools
module load bwa/0.7.7
module load python/2.7.6

## Set directories
    PROJ=/ufrc/mcintyre/share/conesa_isoform_check
    ORIG=$PROJ/original_data
    
    OUTPUT=$PROJ/bwa_mem_aln_pe
    if [ ! -e $OUTPUT ]; then mkdir -p $OUTPUT; fi

    ## Output folder for split SAM files
    SPLITSAM=$OUTPUT/bwa_split
    if [ ! -e $SPLITSAM ]; then mkdir -p $SPLITSAM; fi

    ## Design file
    DESIGN_FILE=$PROJ/design_files/pacbio_samples.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}
    NAME=${SAMPLE}

### (1) Align reads to genome with BWA-MEM

    ## Set reference for BWA MEM
    BWAREF=$PROJ/references/mm10_for_bedtools_v2_BWA

    echo "Aligning reads to genome with BWA-MEM"

    #### run bwa-mem

    bwa mem -t 4 -M $BWAREF $ORIG/${NAME}_SpliceMap.fq > $OUTPUT/${NAME}.sam

### (2) Parse BWA output

    echo "Parsing BWA-MEM SAM files"

    python $PROJ/scripts/BWASplitSAM_4mai.py \
           -fq1 $ORIG/${NAME}_SpliceMap.fq \
            -s $OUTPUT/${NAME}.sam \
            --outdir $SPLITSAM
    ## Cat mapped and opposite
    cat $SPLITSAM/${NAME}_mapped.sam $SPLITSAM/${NAME}_oposite.sam > $OUTPUT/${NAME}_all_uniq.sam

### (3) Convert SAM to BAM, sort and index

    BAMDIR=$PROJ/miso_bam_files_${NAME}
    if [ ! -e $BAMDIR ]; then mkdir -p $BAMDIR; fi

    ## Fasta reference
    FASTA=$PROJ/references/mm10_for_bedtools_v2.fa

    echo "Converting SAM to BAM"

    sam_to_bam --convert $OUTPUT/${NAME}_all_uniq.sam ${BAMDIR} --ref ${FASTA}.fai

### (4) Index GFF file

     echo "Indexing GFF file"

#     GFF=$PROJ/references/aconesa_refseq_gff3_v4.gff
#     GFFDIR=$PROJ/references/gff_indexed
     GFF=$PROJ/references/miso_event_gff/SE.mm10.gff3
     GFFDIR=$PROJ/references/miso_event_gff/gff_indexed_SE
     if [ ! -e $GFFDIR ]; then mkdir -p $GFFDIR; fi

     index_gff --index ${GFF} ${GFFDIR}

### (5) Run MISO

     echo "Running MISO"

     MISO_OUTDIR=$PROJ/miso_output_SE_${NAME}
     if [ ! -e ${MISO_OUTDIR} ]; then mkdir -p ${MISO_OUTDIR}; fi

miso --run ${GFFDIR} $BAMDIR/${NAME}_all_uniq.sorted.bam --output-dir ${MISO_OUTDIR} --read-len 56

     echo "MISO finished"
### (6) Summarize MISO output

     echo "Summarizing MISO output"

     SUMOUT=$PROJ/miso_summary_SE_${NAME}
     if [ ! -e ${SUMOUT} ]; then mkdir -p ${SUMOUT}; fi

     summarize_miso --summarize-samples ${MISO_OUTDIR} ${SUMOUT}
     echo "all done!"
