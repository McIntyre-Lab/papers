#!/bin/bash

#   1. get_components_from_geneID_03mdg.py
#   2. plot_6panel_gene_data.py

GENES=("FBgn0024321")
SOURCE_SPECIES="dmel6"

# PLOT SETTINGS
PANELS="ERPnovel"
OUTPUT_FORMAT="pdf"

ANNO_MIN_READS=1
ERP_ANNO_MIN_READS=10
ERP_ANNO_TOP_N=50
ERP_NOVEL_MIN_READS=50
ERP_NOVEL_TOP_N=50

# PATHS
PROJ="/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO="${PROJ}/zenodo"
DATA="${ZENODO}/datafiles"
ANNO="${ZENODO}/fiveSpecies_supporting_files"
SCRIPTS="${PROJ}/scripts"

COMPONENT_FILE="${ZENODO}/FiveSpecies_network_files/component_map_by_node.csv"
OUTPUT_DIR="${PROJ}/Figures"

GET_COMPONENTS_SCRIPT="${SCRIPTS}/get_components_from_geneID_03mdg.py"
PLOT_SCRIPT="${SCRIPTS}/6panel_plotting_scripts/plot_6panel_gene_data.py"

ANNO_FILE="${ZENODO}/fiveSpecies_${SOURCE_SPECIES}_full_annotation.csv"

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
# PREPARE ARGUMENTS
# ================================
DATA_GTF_ARG="--data_gtf_files ${DATA_GTF_FILES[*]}"
ERP_DATAFILES_ARG="--erp_datafiles ${ERP_DATAFILES[*]}"
ERP_ANNO_TOP_N_ARG="--erp_anno_top_n ${ERP_ANNO_TOP_N}"
ERP_NOVEL_TOP_N_ARG="--erp_novel_top_n ${ERP_NOVEL_TOP_N}"

# ================================
# EXECUTION LOOP
# ================================
mkdir -p "${OUTPUT_DIR}"

for GENE_ID in "${GENES[@]}"; do
    echo "------------------------------------------------"
    echo "Processing Gene: ${GENE_ID}"
    echo "------------------------------------------------"

    # Use gene-specific temp file to avoid race conditions
    TEMP_COMPONENT_LIST="./roz_component_list_${GENE_ID}_$$.csv"

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
        continue 
    fi

    # Step 2: Plot 6-panel figure
    echo "Step 2: Plotting ${GENE_ID} (panels: ${PANELS})..."
    python3 "${PLOT_SCRIPT}" \
        --component_file     "${COMPONENT_FILE}" \
        --component_list     "${TEMP_COMPONENT_LIST}" \
        --gene_name          "${GENE_ID}" \
        --output_dir         "${OUTPUT_DIR}" \
        --panels             "${PANELS}" \
        --output_format      "${OUTPUT_FORMAT}" \
        --anno_min_reads     "${ANNO_MIN_READS}" \
        --erp_anno_min_reads "${ERP_ANNO_MIN_READS}" \
        --erp_novel_min_reads "${ERP_NOVEL_MIN_READS}" \
        --anno_gtf_files     "${ANNO_GTF_FILES[@]}" \
        --jxn_datafiles      "${JXN_DATAFILES[@]}" \
        ${ERP_ANNO_TOP_N_ARG} \
        ${ERP_NOVEL_TOP_N_ARG} \
        ${ERP_ANNO_SIG_ONLY} \
        ${ERP_NOVEL_SIG_ONLY} \
        ${DATA_GTF_ARG} \
        ${ERP_DATAFILES_ARG}

    rm -f "${TEMP_COMPONENT_LIST}"
    echo "Done: ${GENE_ID}"
done

echo "All genes in list have been processed."