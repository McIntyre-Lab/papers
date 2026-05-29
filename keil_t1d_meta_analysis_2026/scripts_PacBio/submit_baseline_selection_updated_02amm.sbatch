#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=base_select
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/baseline_selection_%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --array=2

## Identify and remove unspliced fragments and
##     reads classified as genic
## 1) Flag monoexon reads in classification file
## 2) Flag reads in putative novel/fusion loci
## 3) Flag genic reads
## 4) List reads that are not unspliced fragments (KEEP NOVEL/FUSION AND GENIC)
## 5) Get GTF of isoforms that pass selection parameters
## 6) Subset other outputs for selected reads

module load python/3.6
module load samtools/1.10

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/baseline_selection_keep_genic
    mkdir -p ${OUTD}
IND=$PROJ/sqanti3_qc_allRuns
INDSAM=$PROJ/minimap_w_ref_bed
OUTSAM=$PROJ/minimap_allRuns
    mkdir -p $OUTSAM

## Get info from design file
DESIGN_FILE=$PROJ/design_files/design_sampleID_2_barcode.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
BC=${ARRAY[1]}

SAMPLE_NAME=${SAMPLEID}

date

for RUN in first
do

if [[ ! -e ${OUTSAM}/${SAMPLE_NAME}_minimap2.sam ]]; then
    samtools view -h ${INDSAM}/${SAMPLE_NAME}_minimap2*.bam \
    > ${OUTSAM}/${SAMPLE_NAME}_${RUN}_minimap2.sam
    cp ${INDSAM}/${SAMPLE_NAME}_minimap2*.bam ${OUTSAM}/${SAMPLE_NAME}_${RUN}_minimap2.bam
fi 

done

for RUN in second third fourth fifth
do

if [[ ! -e ${OUTSAM}/${RUN}_run/${SAMPLE_NAME}_minimap2.sam ]]; then
    samtools view -h ${INDSAM}/${RUN}_run/${SAMPLE_NAME}_minimap2*.bam \
    > ${OUTSAM}/${SAMPLE_NAME}_${RUN}_minimap2.sam
    cp ${INDSAM}/${SAMPLE_NAME}_minimap2*.bam ${OUTSAM}/${SAMPLE_NAME}_${RUN}_minimap2.bam

fi
done


for RUN in first second third fourth fifth
do

## Set log file
LOG=${OUTD}/${SAMPLE_NAME}_${RUN}_flag_monoexon_freq.log

echo "
*** Baseline Selection ***

LOG=${OUTD}/${SAMPLE_NAME}_${RUN}_flag_monoexon_freq.log

Inputs: ${IND}/${SAMPLE_NAME}_${RUN}_classification.txt
        ${IND}/${SAMPLE_NAME}_${RUN}_corrected_associated_gene.gtf
"
## Flag monoexon reads in classification file
python ${SCRIPTS}/flag_monoexon_02avn.py \
    -c ${IND}/${SAMPLE_NAME}_${RUN}_classification.txt \
    -o ${OUTD}/${SAMPLE_NAME}_${RUN}_corrected_associated_gene_monoexon_freq.txt \
    > ${LOG}
echo -e "\t...${SAMPLE_NAME} ${RUN} monoexon reads flagged"

## Flag reads in novel/fusion loci or classified as "genic"
## Output list of reads that pass selection parameters
## Keep novel and fusion loci, remove unspliced fragments and genic classified reads
python ${SCRIPTS}/flag_subset_fragment_genic_novel_fusion_reads_03avn.py \
    -c ${IND}/${SAMPLE_NAME}_${RUN}_classification_flag_mono.txt \
    -n \
    -f \
    -g \
    -o ${OUTD}/${SAMPLE_NAME}_${RUN}_corrected_associated_gene_read_flag_freq.txt \
    -p ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt

##rm ${IND}/${SAMPLE_NAME}_${RUN}_classification_flag_mono.txt
echo -e "\t...${SAMPLE_NAME} ${RUN} reads in novel/fusion loci and low evidence genes flagged"

## Get GTF of read that pass selection parameters
python ${SCRIPTS}/subset_gtf_simple_03avn.py \
    -g ${IND}/${SAMPLE_NAME}_${RUN}_corrected_associated_gene.gtf \
    -t transcript_id \
    -i ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt \
    -o ${OUTD}/${SAMPLE_NAME}_${RUN}_corrected_associated_gene_baseline_selection_keep_genic.gtf \
    >> ${LOG}
echo -e "\t...${SAMPLE_NAME} ${RUN} baseline selected GTF output to ${SAMPLE_NAME}_${RUN}_corrected_associated_gene_baseline_selection_keep_genic.gtf"

    ## Count
    echo "
$(awk 'NR!=1' ${IND}/${SAMPLE_NAME}_${RUN}_classification_flag_monoexon_novel_fusion_loci_genic_reads.txt | wc -l) total reads
$(wc -l ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt | awk '{print $1}') pass selection parameters (monoexon reads in multiexon genes or genic classified reads are removed)
">> ${LOG}

## Subset other output files for selected reads
python $PROJ/scripts/Common_extracFeatureFromList.py \
    --list ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt \
    --type class \
    --input ${IND}/${SAMPLE_NAME}_${RUN}_classification_flag_monoexon_novel_fusion_loci_genic_reads.txt \
    --out ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selec_keep_genic_classification.txt

python $PROJ/scripts/Common_extracFeatureFromList.py \
    --list ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt \
    --type fasta \
    --input ${IND}/${SAMPLE_NAME}_${RUN}_corrected.fasta \
    --out ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_select_keep_genic.fasta

python $PROJ/scripts/Common_extracFeatureFromList.py \
    --list ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_selection_keep_genic.txt \
    --type sam \
    --input ${OUTSAM}/${SAMPLE_NAME}_${RUN}_minimap2.sam \
    --out ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_select_keep_genic.sam

samtools view -bS \
    ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_select_keep_genic.sam \
    > ${OUTD}/${SAMPLE_NAME}_${RUN}_baseline_select_keep_genic.bam
#rm ${INDSAM}/${SAMPLE_NAME}_${RUN}_minimap2.sam

done
