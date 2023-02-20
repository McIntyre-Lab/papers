#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --job-name=cvr_frags
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/mcintyre/share/etoh_srna/scripts/SLURM_LOGS/cvr_rna_features_%A-%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1              # single cpu core
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb    		# per processor memory
#SBATCH --array=1-48
#SBATCH --dependency=afterok:41690523
date;hostname;pwd


## Using looping of the array task runs
## df (no header) has 48 rows
## want 1 rows per array task
        ## split into X array tasks:  48/X = 1 rows per array task, X=48
#Set the number of runs that each SLURM task should do
PER_TASK=1

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo "This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM"

# Run the loop of runs for this task.
for (( RUN=$START_NUM; RUN<=$END_NUM; RUN++ ))
do

    echo "This is SLURM task $SLURM_ARRAY_TASK_ID, run number $RUN"

    module purge
    module load gcc/5.2.0
    module load python/2.7.10
    module load samtools/1.9

    ## Set directories
    PROJ=/ufrc/mcintyre/share/etoh_srna
    SCRIPTS=/ufrc/mcintyre/share/python.git
    ALN=$PROJ/rnaseq/aln_genome
    BAMFILES=${ALN}/combine_techRep_bam_files
    PILES=${ALN}/mpileups
        mkdir -p ${PILES}

    ## Design file                       
    DESIGN_FILE=$PROJ/design_files/design_ethanol_srna_RNAseq_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SPECIES=${ARRAY[1]}
    GENO=${ARRAY[2]}
    SEX=${ARRAY[3]}
    TRT=${ARRAY[4]}
    REP=${ARRAY[5]}
    WELLID=${ARRAY[7]}

    SAMPLEID=${SPECIES}_${GENO}_${SEX}_${TRT}_rep${REP}
    echo "${SAMPLEID}...
"
    if [[ ${SPECIES} == "mel" ]]; then
        REF=/ufrc/mcintyre/share/references/dmel_fb617/dmel-all-chromosome-r6.17.fasta
        FRAGMENTS=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_exon_fragments_coverage.bed
        FUSIONS=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_fusions_coverage.bed
        INTRONS=/ufrc/mcintyre/share/references/dmel_fb617/dmel617_introns_from_fusions.bed
        echo "Species is ${SPECIES} and ref is $REF"
    else
        REF=/ufrc/mcintyre/share/references/dsim_fb202/dsim-all-chromosome-r2.02.fasta
        FRAGMENTS=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_150bp_reads/dsim202_150bp_exon_fragments_coverage.bed
        FUSIONS=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_150bp_reads/dsim202_150bp_fusions_coverage.bed
        INTRONS=/ufrc/mcintyre/share/references/dsim_fb202/dsim202_annotations_150bp_reads/dsim202_150bp_introns_from_fusions.bed
        echo "Species is ${SPECIES} and ref is ${REF}"
    fi
    TSS=$PROJ/dros_annotation_pipeline/${SPECIES}_TSS300bpWindow_unique.bed
    fiveUTR=$PROJ/dros_annotation_pipeline/${SPECIES}_5UTR_unique.bed
    threeUTR=$PROJ/dros_annotation_pipeline/${SPECIES}_3UTR_unique.bed
    INTERGENIC=$PROJ/dros_annotation_pipeline/${SPECIES}_intergenic.bed


    ## Create mpileup from merged bam files
    if [ ! -e ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup ]; then
        samtools mpileup \
            -d 1000000 \
            -f $REF \
            ${BAMFILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam \
            > ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup
    fi

    ## convert sorted bam to sam
    if [ ! -e ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam ]; then
        samtools view -h -o ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam ${BAMFILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.bam
    fi

    for BED in ${FRAGMENTS} ${FUSIONS} ${INTRONS} ${TSS} ${fiveUTR} ${threeUTR} ${INTERGENIC}; do

        if [[ ${BED} == ${FRAGMENTS} ]]; then
            TYPE="fragment"
        elif [[ ${BED} == ${FUSIONS} ]]; then
            TYPE="fusion"
        elif [[ ${BED} == ${INTRONS} ]]; then
            TYPE="intron"
        elif [[ ${BED} == ${TSS} ]]; then
            TYPE="TSS300bpWindow"
        elif [[ ${BED} == ${fiveUTR} ]]; then
            TYPE="5UTR"
        elif [[ ${BED} == ${threeUTR} ]]; then
            TYPE="3UTR"
        elif [[ ${BED} == ${INTERGENIC} ]]; then
            TYPE="intergenic"
        else
            echo "ERROR: UNKNOWN BED FILE"
            exit 1
        fi

        ## Set output directories
        CVR=$PROJ/rnaseq/coverage_cnts/${TYPE}_coverage_cnts
            mkdir -p ${CVR}
        LOGC=${CVR}/logs
            mkdir -p ${LOGC}

        ## Get coverage counts if they do not already exist
        echo -e "\t${TYPE} coverage counts..."
        if [ ! -e ${CVR}/${feature}/cvrg_cnts_${SAMPLEID}_combined.csv ]; then
            ## get coverage counts for alignments using rpkm_calculate
            python $SCRIPTS/rpkm_calculate.py \
                -m ${PILES}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.mpileup \
                -s ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam \
                -b ${BED} \
                -g ${LOGC}/cov_cnts_aln_${SAMPLEID}_combined.logfile \
                -o ${CVR}/cvrg_cnts_${SAMPLEID}_combined.csv \
                -c
        fi
    done

    ## remove sam and mpileup files
    rm ${ALN}/${SAMPLEID}_combined_bwa_PE_SE_uniq.sorted.sam

done
