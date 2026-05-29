#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=trand_isoseq.%j
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/trand_isoseq.%j
#SBATCH --cpus-per-task=64
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-25


###Set directories and files
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
SQANTI=${PROJ}/sqanti3_isoseq_filtered_transcripts_vs_union_02ndk
OUTD=${PROJ}/trand_gene_isoseq_FSM_ISM_NIC_02ndk
    mkdir -p ${OUTD}
ROZ=$OUTD/roz
    mkdir -p ${ROZ}

## design does not have header!!
DESIGN_FILE=$PROJ/design_files/design_chromosomes.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

CHROM=${ARRAY[0]}


REF_GTF=/blue/mcintyre/share/references/hg38_union/hg38_union_UJC.gtf

GTF=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk/isoseq_allSamples_corrected_geneID.gtf

module load python/3.6

####Select reference transcripts to FSMs and ISMs and select NIC transcript models


##Selecting reference transcripts for FSMs and ISMs to union ref

   awk '$6 == "full-splice_match" || $6 == "incomplete-splice_match" {print $8}' \
        ${SQANTI}/${CHROM}_isoseq_vs_union_classification.txt | \
        sort | uniq \
        >${ROZ}/${CHROM}_FSM_ISM_reference_transcripts.txt

    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${REF_GTF} \
        -t transcript_id \
        -i ${ROZ}/${CHROM}_FSM_ISM_reference_transcripts.txt \
        -o ${ROZ}/${CHROM}_FSM_ISM_reference_transcripts.gtf



##Selecting isoseq transcript models that are NICs to union ref

    awk '$6 == "novel_in_catalog" {print $1}' \
        ${SQANTI}/${CHROM}_isoseq_vs_union_classification.txt \
        >${ROZ}/${CHROM}_NIC_isoseq_transcript_models.txt

    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${GTF} \
        -t transcript_id \
        -i ${ROZ}/${CHROM}_NIC_isoseq_transcript_models.txt \
        -o ${ROZ}/${CHROM}_NIC_isoseq_transcript_models.gtf


cat ${ROZ}/${CHROM}_FSM_ISM_reference_transcripts.gtf ${ROZ}/${CHROM}_NIC_isoseq_transcript_models.gtf > ${OUTD}/${CHROM}_isoseq_FSM_ISM_NIC.gtf

module purge
export PATH=/blue/mcintyre/ammorse/conda_envs/trand_dev/bin:$PATH

date
echo "
*** TranD pairwise isoseq transcripts ***
"
    trand \
	${OUTD}/${CHROM}_isoseq_FSM_ISM_NIC.gtf \
        -e gene \
        -k \
	-l ${OUTD}/${CHROM}_isoseq_FSM_ISM_NIC.log \
        -o ${OUTD} \
        -x ${CHROM}_isoseq_FSM_ISM_NIC \
        -f \
	-n $SLURM_CPUS_ON_NODE
    date

