#!/bin/bash

# ================================
# INPUT
# ================================
#GENESETID="23"
COMPONENT_IDS=(
    7295
    36119
)
# ================================
# PLOT SETTINGS
# ================================
PANELS="ERPnovel"
OUTPUT_FORMAT="pdf"

ANNO_MIN_READS=0
ERP_ANNO_MIN_READS=10
ERP_ANNO_TOP_N=""
ERP_NOVEL_MIN_READS=50
ERP_NOVEL_TOP_N=""

# Annotated UJCs with 0 reads are retained and drawn as gray outlines.
# Mono-exon UJCs are dropped by default.

# ================================
# PATHS
# ================================
PROJ="/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO="${PROJ}/zenodo"
DATA="${ZENODO}/datafiles"
ANNO="${ZENODO}/fiveSpecies_supporting_files"
SCRIPTS="${PROJ}/scripts"
OUTPUT_DIR="${PROJ}/Figures"

PLOT_SCRIPT="${SCRIPTS}/6panel_plotting_scripts/plot_6panel_gene_data.py"

FULL_ANNO_FILES=(
    "dmel6:${ZENODO}/fiveSpecies_dmel6_full_annotation.csv"
    "dsim2:${ZENODO}/fiveSpecies_dsim2_full_annotation.csv"
    "dyak2:${ZENODO}/fiveSpecies_dyak2_full_annotation.csv"
    "dsan1:${ZENODO}/fiveSpecies_dsan1_full_annotation.csv"
    "dser1:${ZENODO}/fiveSpecies_dser1_full_annotation.csv"
)

ANNO_GTF_FILES=(
    "dmel6:${ANNO}/fiveSpecies_2_dmel6_anno_files/fiveSpecies_2_dmel6_ujc.gtf"
    "dsim2:${ANNO}/fiveSpecies_2_dsim2_anno_files/fiveSpecies_2_dsim2_ujc.gtf"
    "dyak2:${ANNO}/fiveSpecies_2_dyak2_anno_files/fiveSpecies_2_dyak2_ujc.gtf"
    "dsan1:${ANNO}/fiveSpecies_2_dsan1_anno_files/fiveSpecies_2_dsan1_ujc.gtf"
    "dser1:${ANNO}/fiveSpecies_2_dser1_anno_files/fiveSpecies_2_dser1_ujc.gtf"
)

JXN_DATAFILES=(
    "dmel6:${DATA}/datafile_jxnHash_dmel6.csv"
    "dsim2:${DATA}/datafile_jxnHash_dsim2.csv"
    "dyak2:${DATA}/datafile_jxnHash_dyak2.csv"
    "dsan1:${DATA}/datafile_jxnHash_dsan1.csv"
    "dser1:${DATA}/datafile_jxnHash_dser1.csv"
)

DATA_GTF_FILES=(
    "dmel6:${DATA}/dmel_data_2_dmel6_ujc_noMultiGene.gtf"
    "dsim2:${DATA}/dsim_data_2_dsim2_ujc_noMultiGene.gtf"
    "dyak2:${DATA}/dyak_data_2_dyak2_ujc_noMultiGene.gtf"
    "dsan1:${DATA}/dsan_data_2_dsan1_ujc_noMultiGene.gtf"
    "dser1:${DATA}/dser_data_2_dser1_ujc_noMultiGene.gtf"
)

ERP_DATAFILES=(
    "dmel6:${DATA}/datafile_erp_dmel6.csv"
    "dsim2:${DATA}/datafile_erp_dsim2.csv"
    "dyak2:${DATA}/datafile_erp_dyak2.csv"
    "dsan1:${DATA}/datafile_erp_dsan1.csv"
    "dser1:${DATA}/datafile_erp_dser1.csv"
)

# ================================
# EXECUTION
# ================================
mkdir -p "${OUTPUT_DIR}"

ERP_ANNO_TOP_N_ARG=""
if [[ -n "${ERP_ANNO_TOP_N}" ]]; then
    ERP_ANNO_TOP_N_ARG="--erp_anno_top_n ${ERP_ANNO_TOP_N}"
fi

ERP_NOVEL_TOP_N_ARG=""
if [[ -n "${ERP_NOVEL_TOP_N}" ]]; then
    ERP_NOVEL_TOP_N_ARG="--erp_novel_top_n ${ERP_NOVEL_TOP_N}"
fi

echo "------------------------------------------------"
echo "Processing components: ${COMPONENT_IDS[*]}"
echo "------------------------------------------------"

python3 "${PLOT_SCRIPT}" \
    --component_ids "${COMPONENT_IDS[@]}" \
    --full_anno_files     "${FULL_ANNO_FILES[@]}" \
    --output_dir          "${OUTPUT_DIR}" \
    --panels              "${PANELS}" \
    --output_format       "${OUTPUT_FORMAT}" \
    --anno_min_reads      "${ANNO_MIN_READS}" \
    --erp_anno_min_reads  "${ERP_ANNO_MIN_READS}" \
    --erp_novel_min_reads "${ERP_NOVEL_MIN_READS}" \
    ${ERP_ANNO_TOP_N_ARG} \
    ${ERP_NOVEL_TOP_N_ARG} \
    --anno_gtf_files      "${ANNO_GTF_FILES[@]}" \
    --jxn_datafiles       "${JXN_DATAFILES[@]}" \
    --data_gtf_files      "${DATA_GTF_FILES[@]}" \
    --erp_datafiles       "${ERP_DATAFILES[@]}"

echo "Done: components ${COMPONENT_IDS[*]}"
