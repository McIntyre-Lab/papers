#!/bin/sh
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --job-name=frag_cvr5
#SBATCH --mail-type=FAIL,END
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/frag_cvr5.%A-%a.out
#SBATCH --time=72:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb    		# per processor memory
#SBATCH --array=2-56

date;hostname;pwd

### create mpileups and fragment coverage counts for reads aligned to hg38  genome

module load bwa/0.7.7 gcc/5.2.0 python/2.7.10 samtools/1.12

## Set directories
PROJ=/blue/concannon/share/T1DGC_R01_controls_RNAseq
SCRIPTS=$PROJ/scripts
ALN=$PROJ/alignments_hg38
ALN_COMB=${ALN}/combine_techReps
    mkdir -p ${ALN_COMB}
BAMFILES=$ALN

PILES=$ALN_COMB/mpileups
    mkdir -p $PILES
CVR=$PROJ/coverage_cnts_FSM_ISM_NIC_frags
    mkdir -p $CVR
LOGC=$CVR/logs
    mkdir -p $LOGC

## Design file

    DESIGN_FILE=$PROJ/design_files/design_set5_sampleID.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"


    SAMPLE=${ARRAY[0]}

    REF=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa
    BED=/blue/concannon/share/T1D_R01_controls_PacBio/trand_gene_isoseq_FSM_ISM_NIC_02ndk/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.bed

   ##combine bam files

  samtools merge $ALN_COMB/${SAMPLE}.bam $BAMFILES/${SAMPLE}_*_sorted.bam 


   ## sort bam files
   if [[  ! -e $ALN_COMB/${SAMPLE}_sorted.bam ]]
   then
   samtools sort -T ${SAMPLE}.tmp.sorted.bam -o $ALN_COMB/${SAMPLE}_sorted.bam $ALN_COMB/${SAMPLE}.bam
   samtools index $ALN_COMB/${SAMPLE}_sorted.bam
   fi

   ## create mpileup from bam files
    samtools mpileup \
        -d 1000000 \
        -f $REF \
        $ALN_COMB/${SAMPLE}_sorted.bam \
        > $PILES/${SAMPLE}_sorted.mpileup

    ## convert sorted bam to sam
    samtools view -h -o $ALN_COMB/${SAMPLE}_sorted.sam $ALN_COMB/${SAMPLE}_sorted.bam

    ## get coverage counts for alignments using rpkm_calculate
    python $SCRIPTS/rpkm_calculate.py \
        -m $PILES/${SAMPLE}_sorted.mpileup \
        -s $ALN_COMB/${SAMPLE}_sorted.sam \
        -b $BED \
        -g $LOGC/cvrg_${SAMPLE}.logfile \
        -o $CVR/cvrg_${SAMPLE}.csv \
        -c

    rm $PILES/${SAMPLE}_sorted.mpileup
    rm $ALN_COMB/${SAMPLE}_sorted.sam

    if [[ -e $ALN_COMB/${SAMPLE}.bam ]]
    then
    rm $ALN_COMB/${SAMPLE}.bam
    fi
