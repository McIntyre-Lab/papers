#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=plot
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/plot_transcriptome_MaizeWang2018_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --array=6

module load python/3.7

## Plot transcripts per gene, number of exon regions per gene,
##     number of exonic fragments per gene, and proportion of
##     varying exonic regions/fragments per gene for consolidated
##     representative reference transcripts of FSM/ISM
##     and the monoexon filtered NIC/NNC

## Get info from design file
PROJ=/blue/mcintyre/share/maize_ainsworth
DESIGN_FILE=$PROJ/design_files/df_MaizeWang2018_maize_pool1_sample_2_barcode_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

BC=${ARRAY[0]}
SPECIES=${ARRAY[1]}
TISSUE=${ARRAY[2]}

## Get sampleID (species_tissue)
SAMPLEID=${SPECIES}_${TISSUE}

## Set directories
SCRIPTS=$PROJ/scripts/pacbio_transcriptome_eval
DIST=$PROJ/MaizeWang2018_transcriptome_eval/transcriptDistance/${SAMPLEID}_transcriptDistance
IND=$PROJ/MaizeWang2018_transcriptome_eval/EA_output_curated_transcriptome/${SAMPLEID}_150bp_annotations
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome
    mkdir -p ${OUTD}

## Get list of transcripts not present in EA annotations
awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
    ${IND}/${SAMPLEID}_event2transcript2gene_index.csv | \
    sed s/'|'/'\n'/g | sort | uniq > ${OUTD}/${SAMPLEID}_ea_output_transcript.txt
awk '$3=="exon"{if($9=="transcript_id"){print $10;}if($11=="transcript_id"){print $12}}' \
    $PROJ/MaizeWang2018_transcriptome_eval/subset_FSM_ISM_NIC_NNC/${SAMPLEID}_MaizeWang2018_fsm_ism_nic_nnc.gtf \
    | sed s/'\"'/''/g | sed s/';'/''/g | sort | uniq > ${OUTD}/${SAMPLEID}_input_transcript.txt
diff ${OUTD}/${SAMPLEID}_input_transcript.txt ${OUTD}/${SAMPLEID}_ea_output_transcript.txt |\
    grep "<" |  sed s/"< "/""/g > ${OUTD}/${SAMPLEID}_diff_input_ea_output_transcript.csv
#rm ${OUTD}/${SAMPLEID}_input_transcript.txt ${OUTD}/${SAMPLEID}_ea_output_transcript.txt

## Subset EA annotations for consolidated transcriptome
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/${SAMPLEID}_event2transcript2gene_index.csv \
#    -v event_id \
#    -l ${DIST}/${SAMPLEID}_transcript_list_not_retained.csv \
#    -o ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/${SAMPLEID}_fusion_annotations.csv \
#    -v fusion_id \
#    -l ${DIST}/${SAMPLEID}_transcript_list_not_retained.csv \
#    -o ${IND}/${SAMPLEID}_consol_fusion_annotations.csv
#python ${SCRIPTS}/subset_EA_annotation_files.py \
#    -i ${IND}/${SAMPLEID}_exon_fragment_annotations.csv \
#    -v fragment_id \
#    -l ${DIST}/${SAMPLEID}_transcript_list_not_retained.csv \
#    -o ${IND}/${SAMPLEID}_consol_exon_fragment_annotations.csv

## Print total transcripts and genes represented
##     NIC/NNC transcripts have "PB" in the name,
##     FSM/ISM representative reference transcripts do not
#echo "
#$(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
#    ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv | \
#    sed s/'|'/'\n'/g | sort | uniq | wc -l) transcripts represented ($(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
#    ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv | \
#    sed s/'|'/'\n'/g | sort | uniq | grep -v "PB" | wc -l) consolidated FSM/ISM representative reference transcripts and $(awk -F "," 'NR!=1&&$2!="Unannotated"{print $2}' \
#    ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv | \
#    sed s/'|'/'\n'/g | sort | uniq | grep "PB" | wc -l) consolidated filtered NIC/NNC)
#$(awk -F "," 'NR!=1{print $3}' \
#    ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv | \
#    sed s/'|'/'\n'/g | sort | uniq | wc -l) reference genes represented
#" > ${OUTD}/${SAMPLEID}_consol_transcript_gene_counts.txt

## Plot
#python ${SCRIPTS}/EA_histgram_visualize_07lzh.py \
#    -ig ${IND}/${SAMPLEID}_consol_event2transcript2gene_index.csv \
#    -ie ${IND}/${SAMPLEID}_consol_fusion_annotations.csv \
#    -if ${IND}/${SAMPLEID}_consol_exon_fragment_annotations.csv \
#    --prefix ${SAMPLEID}_consol \
#    -o ${OUTD} \
#    -s -l -r
