#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=sqanti_isoseq.%j
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sqanti_isoseq.%j
#SBATCH --cpus-per-task=64
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-25

module purge
module load python/3.6
module load sqanti3/5.1.1

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
IND=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk
OUTD=${PROJ}/sqanti3_isoseq_filtered_transcripts_vs_union_02ndk
    mkdir -p ${OUTD}
ROZ=$OUTD/roz
    mkdir -p ${ROZ}


## design does not have header!!
DESIGN_FILE=$PROJ/design_files/design_chromosomes.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

CHROM=${ARRAY[0]}


SQANTI=${PROJ}/sqanti3_qc_isoseq_vs_UJCs
MONO=${PROJ}/flag_monoexon_isoseq_allSamples_02ndk/flag_monoexon_isoseq_allSamples.csv
GTF=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk/isoseq_allSamples_corrected_geneID.gtf

echo "
*** Selecting isoseq transcripts that are FSMs to the consolidated reads**
"
    ## Get isoseq transcripts that are FSMs to the consolidated reads
    awk '$6 == "full-splice_match" {print $1}' \
        ${SQANTI}/${CHROM}_classification.txt \
        >${ROZ}/${CHROM}_list_FSM_to_reads.txt

echo "
*** Selecting isoseq transcripts that are unsspliced fragments**
"

    awk -F "," '$7 == 1 {print $4}' \
        ${MONO} \
        >${ROZ}/list_unspliced_fragments.txt

echo "
*** Subset isoseq gtf for reads that are FSMs to reads AND are NOT unspliced fragments**
"

    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${GTF} \
        -t transcript_id \
        -i ${ROZ}/${CHROM}_list_FSM_to_reads.txt \
        -o ${ROZ}/${CHROM}_subset_FSM_2_reads.gtf


    python $PROJ/scripts/subset_gtf_simple_03avn.py \
        -g ${ROZ}/${CHROM}_subset_FSM_2_reads.gtf \
        -t transcript_id \
        -e ${ROZ}/list_unspliced_fragments.txt \
        -o ${OUTD}/${CHROM}_filtered_isoseq_trancript_models.gtf



#After removing transcript models that are not FSMs to the reads and unspliced fragments
#Run sqanti3 qc of these transcript models against the union reference

REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa
REF_GTF=/blue/mcintyre/share/references/hg38_union/hg38_union_UJC.gtf


##  --gtf option indicates GTF isoform input
##  --skipORF option to skip ORF prediction (runs much faster)
##  --force_id_ignore option to allow transcript ID names that are not like PB.X.Y
##  Order of files is (1) isoforms (FA/FQ/GTF), (2) annotation (ref GTF/GFF), (3) genome (ref FA)

#   --dir is path to output directory
#   -o prefix name to use

## order of files input into sqanti3_qc.py is critical!!
# 1. sam2gtf (from read align to genome)
# 2. ref gtf file
# 3. ref genome fasta file


sqanti3_qc.py \
    --skipORF \
    --force_id_ignore \
    --dir ${OUTD} \
    --report skip \
    -o ${CHROM}_isoseq_vs_union \
    ${OUTD}/${CHROM}_filtered_isoseq_trancript_models.gtf \
    ${REF_GTF} \
    ${REF_FA}



