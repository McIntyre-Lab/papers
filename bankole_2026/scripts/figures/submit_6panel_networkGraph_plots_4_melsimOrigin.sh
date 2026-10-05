#!/bin/bash

# Plot network graphs for selected five-species components.

###############################################################################
# directory paths

PROJ="/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing"
ZENODO="${PROJ}/zenodo"
SCRIPTS="${PROJ}/scripts"
NETWORK_SCRIPT="${SCRIPTS}/plot_component_network_graph_01mdg.py"

INPUT_DIR="${ZENODO}/FiveSpecies_network_files"
OUTPUT_DIR="${PROJ}/Figures"

###############################################################################
# inputs & arguments

COMPONENT_IDS=(
#    7632
#    14324
#    14729
#    4282
#    15471
#    11478
#    14049
#    14154
#    9327
#    5762
#    14785
#    15639
#    2774
#    7024
#    8992
#    1089
#    10265
#    3647
#    14151
#    718
#    14025
#    15178
#    15756
#    9568
#    8478
#    9685
#    2424
#    332
    14886
)

NODES_FILE="${INPUT_DIR}/component_map_by_node.csv"
EDGES_FILE="${INPUT_DIR}/edges.csv"

# ================================
# EXECUTION
# ================================

mkdir -p "${OUTPUT_DIR}"

for COMPONENT_ID in "${COMPONENT_IDS[@]}"
do

    echo "------------------------------------------------"
    echo "Processing component ${COMPONENT_ID}"
    echo "------------------------------------------------"

    python3 "${NETWORK_SCRIPT}" \
        --component_id "${COMPONENT_ID}" \
        --nodes_file   "${NODES_FILE}" \
        --edges_file   "${EDGES_FILE}" \
        --output_dir   "${OUTPUT_DIR}"

    echo "Done: component ${COMPONENT_ID}"

done
