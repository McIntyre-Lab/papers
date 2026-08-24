#!/bin/bash
set -euo pipefail

PROJ=/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing
SCRIPTS=${PROJ}/scripts
OUT=${PROJ}/zenodo/summary_files
ROZ=/TB14/TB14/mgaran/roz

GROUP_CLASS_COL="sexClass"
CPUS=24

echo "=========================================="
echo "AS event detection (standalone)"
echo "Group class column: $GROUP_CLASS_COL"
echo "CPUs for AS detection: $CPUS"
echo "Reads flagged UJC/ERP propReads tables left in ROZ by the gene summary pipeline"
echo "=========================================="

DESIGN=${PROJ}/design_files/fiveSpecies_M_F_designFile.csv

while IFS=',' read -r SPECIES GROUP1 GROUP2 COUNT_PREFIX ; do
    SPECIES=$(echo "$SPECIES" | tr -d '\r' | xargs)
    GROUP1=$(echo "$GROUP1" | tr -d '\r' | xargs)
    GROUP2=$(echo "$GROUP2" | tr -d '\r' | xargs)

    if [ -z "$SPECIES" ]; then
        continue
    fi

    echo ""
    echo "=========================================="
    echo "Processing $SPECIES"
    echo "=========================================="

    UJC_FLAGGED="$ROZ/min10reads_UJC_propReads_table_w_flags_${SPECIES}.csv"
    ERP_FLAGGED="$ROZ/min10reads_ERP_propReads_table_w_flags_${SPECIES}.csv"
    AS_GENE_SUM="$OUT/gene_summary_min10reads_from_AS_analysis_${SPECIES}.csv"

    echo "Detecting AS events (standalone table)"
    python "$SCRIPTS/detect_AS_events.py" \
        -i "$UJC_FLAGGED" -e "$ERP_FLAGGED" -o "$AS_GENE_SUM" -s "$SPECIES" \
        -g1 "$GROUP1" -g2 "$GROUP2" --group-class-col "$GROUP_CLASS_COL" --cpus "$CPUS"

    echo "Completed $SPECIES"
done < <(tail -n +2 "$DESIGN")

echo ""
echo "AS event detection complete"
echo "Final outputs in: $OUT"
echo "  gene_summary_min10reads_from_AS_analysis_{species}.csv"
