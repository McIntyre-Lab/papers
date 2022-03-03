#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=htseq_count
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/htseq_count_genes_loci_%A-%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --array=1-120

# Using looping of the array task runs
## df (no header) has 120 rows

#Set the number of runs that each SLURM task should do
##PER_TASK=5
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
    module load htseq/0.11.2

    ### Set Directories
    PROJ=/blue/mcintyre/share/maize_ainsworth
    IND=$PROJ/aln_rnaseq_bwa_mem_genome
    OUTD=$PROJ/htseq_quantify_genes_loci/htseq_output_counts
        mkdir -p ${OUTD}

    ## Design file
    DESIGN_FILE=$PROJ/design_files/maize_rnaseq_samples_fix_noHeader.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $RUN | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    SAMPLE=${ARRAY[0]}

    ## GTF for quantification
    GTF=$PROJ/htseq_quantify_genes_loci/sqanti_b73_filtered_corrected_associated_gene.gtf

    ## Quantify genes/loci with htseq
    ## SE
    htseq-count \
        ${IND}/${SAMPLE}_merge_bwa_SE.sam \
        ${GTF} \
        -f sam \
        -s no \
        -m union \
        --secondary-alignment=ignore \
        --supplementary-alignments=ignore \
        > ${OUTD}/${SAMPLE}_SE_htseq_count.tsv

    ## PE
    htseq-count \
        ${IND}/${SAMPLE}_unmerge_bwa_PE.sam \
        ${GTF} \
        -f sam \
        -s no \
        -m union \
        --secondary-alignment=ignore \
        --supplementary-alignments=ignore \
        > ${OUTD}/${SAMPLE}_PE_htseq_count.tsv

done
