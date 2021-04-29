#!/bin/bash
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=plot
#SBATCH	--output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/plot_vs_ref_MaizeWang2018_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-6

module load python/3.7

## Plot transcripts per gene for consolidated
##     representative reference transcripts of FSM/ISM
##     and the monoexon filtered NIC/NNC  vs. reference


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
OUTD=$PROJ/MaizeWang2018_transcriptome_eval/plot_consol_curated_transcriptome_vs_ref
    mkdir -p ${OUTD}

DATANAME="Wang, et al. 2018 ${SAMPLEID}"

## Get FSM consolidated reference EA event2transcript2gene file
if [[ ${SPECIES} == "maize" ]]; then
    REFEVENT=/blue/mcintyre/share/references/maize_b73/FSM_consolidation_maize_B73_EA_150bp/FSM_consolidation_maize_B73_event2transcript2gene_index.csv
    REFNAME="Maize B73 Reference"
else
    REFEVENT=/blue/mcintyre/share/references/sorghum_bicolor/FSM_consolidation_sorghum_EA_150bp/FSM_consolidation_sorghum_event2transcript2gene_index.csv
    REFNAME="Sorghum Bicolor Reference"
fi

date

echo "
*** Plot comparison of data to reference ***
"
EAPREFIX=${IND}/${SAMPLEID}

python ${SCRIPTS}/plot_compare_xcrpt_per_gene_02avn.py \
    -1 ${REFEVENT} \
    -2 ${EAPREFIX}_event2transcript2gene_index.csv \
    -n1 "${REFNAME}" \
    -n2 "${DATANAME}" \
    --hist-1-2 \
    --scatter \
    -d ${OUTD} \
    -p ${SAMPLEID}_vs_reference \
    > ${OUTD}/${SAMPLEID}_compare_2_ref_counts.txt

