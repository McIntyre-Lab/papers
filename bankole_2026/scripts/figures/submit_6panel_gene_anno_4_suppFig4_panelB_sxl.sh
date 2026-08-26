#!/bin/bash
#
# Produces an annotation-only species-color 6-panel plot for a gene
# and its component network graph.
#
# Steps executed:
#   1. get_components_from_geneID_03mdg.py
#   2. plot_6panel_gene_anno.py

# INPUT
GENE_ID="FBgn0264270"
SOURCE_SPECIES="dmel6"

# PLOT SETTINGS
OUTPUT_FORMAT="png"   # svg and pdf are other options

# PATHS
PROJ="/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO="${PROJ}/zenodo"
DATA="${ZENODO}/datafiles"
ANNO="${ZENODO}/fiveSpecies_supporting_files"
SCRIPTS="${PROJ}/scripts"

COMPONENT_FILE="${ZENODO}/fiveSpecies_network_files/component_map_by_node.csv"
OUTPUT_DIR="${PROJ}/Figures/6panel_plots/gene_anno"

GET_COMPONENTS_SCRIPT="${SCRIPTS}/get_components_from_geneID_03mdg.py"
PLOT_SCRIPT="${SCRIPTS}/6panel_plotting_scripts/plot_6panel_gene_anno.py"
NETWORK_SCRIPT="${SCRIPTS}/plot_component_network_graph_01mdg.py"

ANNO_FILE="${ZENODO}/fiveSpecies_${SOURCE_SPECIES}_full_annotation.csv"

ANNO_GTF_FILES=(
    "dmel6:${ANNO}/fiveSpecies_2_dmel6_anno_files/fiveSpecies_2_dmel6_ujc.gtf"
    "dsim2:${ANNO}/fiveSpecies_2_dsim2_anno_files/fiveSpecies_2_dsim2_ujc.gtf"
    "dyak2:${ANNO}/fiveSpecies_2_dyak2_anno_files/fiveSpecies_2_dyak2_ujc.gtf"
    "dsan1:${ANNO}/fiveSpecies_2_dsan1_anno_files/fiveSpecies_2_dsan1_ujc.gtf"
    "dser1:${ANNO}/fiveSpecies_2_dser1_anno_files/fiveSpecies_2_dser1_ujc.gtf"
)

FULL_ANNO_FILES=(
    "dmel6:${ZENODO}/fiveSpecies_dmel6_full_annotation.csv"
    "dsim2:${ZENODO}/fiveSpecies_dsim2_full_annotation.csv"
    "dyak2:${ZENODO}/fiveSpecies_dyak2_full_annotation.csv"
    "dsan1:${ZENODO}/fiveSpecies_dsan1_full_annotation.csv"
    "dser1:${ZENODO}/fiveSpecies_dser1_full_annotation.csv"
)

# ================================
# EXECUTION
# ================================
mkdir -p "${OUTPUT_DIR}"
TEMP_COMPONENT_LIST="./roz_component_list_$$.csv"

# Step 1: Resolve gene to component list
echo "Step 1: Getting components for ${GENE_ID} (${SOURCE_SPECIES})..."
python3 "${GET_COMPONENTS_SCRIPT}" \
    --component_file "${COMPONENT_FILE}" \
    --target_gene    "${GENE_ID}" \
    --source_species "${SOURCE_SPECIES}" \
    --output         "${TEMP_COMPONENT_LIST}" \
    --add_erp        "${SOURCE_SPECIES}" \
    --anno_file      "${ANNO_FILE}"

if [[ ! -f "${TEMP_COMPONENT_LIST}" ]]; then
    echo "ERROR: Failed to generate component list for ${GENE_ID}"
    exit 1
fi

# Step 2: Plot annotation-only species-color 6-panel figure
echo "Step 2: Plotting ${GENE_ID}..."
python3 "${PLOT_SCRIPT}" \
    --component_file  "${COMPONENT_FILE}" \
    --component_list  "${TEMP_COMPONENT_LIST}" \
    --gene_name       "${GENE_ID}" \
    --output_dir      "${OUTPUT_DIR}" \
    --anno_gtf_files  "${ANNO_GTF_FILES[@]}" \
    --full_anno_files "${FULL_ANNO_FILES[@]}" \
    --output_format   "${OUTPUT_FORMAT}"

rm -f "${TEMP_COMPONENT_LIST}"
echo "Done: ${GENE_ID}"