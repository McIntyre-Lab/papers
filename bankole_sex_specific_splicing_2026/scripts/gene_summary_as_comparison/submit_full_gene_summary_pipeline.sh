#!/bin/bash
set -euo pipefail

PROJ=/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing
DATA=${PROJ}/zenodo/datafiles
ANNO=${PROJ}/zenodo/fiveSpecies_supporting_files
SCRIPTS=${PROJ}/scripts
OUT=${PROJ}/zenodo/summary_files
ROZ=/TB14/TB14/mgaran/roz
    mkdir -p ${ROZ}
LOG=${ROZ}/log_files
    mkdir -p ${LOG}

UJC_ANNOTATION_FLAG=flag_jxnHash_in_fiveSpecies_full_anno
ERP_ANNOTATION_FLAG=flag_ERPanno
GROUP_NAME="Sex"
GROUP_CLASS_COL="sexClass"
CPUS=24
MIN_READS=1

declare -A GENOME_MAP=(
    [dmel]="dmel6"
    [dsim]="dsim2"
    [dser]="dser1"
    [dsan]="dsan1"
    [dyak]="dyak2"
)

echo "=========================================="
echo "Full gene summary pipeline"
echo "Group class column: $GROUP_CLASS_COL"
echo "min-reads threshold: $MIN_READS"
echo "CPUs for AS detection: $CPUS"
echo "=========================================="

DESIGN=${PROJ}/design_files/fiveSpecies_M_F_designFile.csv

while IFS=',' read -r SPECIES GROUP1 GROUP2 COUNT_PREFIX ; do
    SPECIES=$(echo "$SPECIES" | tr -d '\r' | xargs)
    GROUP1=$(echo "$GROUP1" | tr -d '\r' | xargs)
    GROUP2=$(echo "$GROUP2" | tr -d '\r' | xargs)
    COUNT_PREFIX=$(echo "$COUNT_PREFIX" | tr -d '\r' | xargs)

    if [ -z "$SPECIES" ]; then
        continue
    fi

    echo ""
    echo "=========================================="
    echo "Processing $SPECIES"
    echo "=========================================="

    GENOME=${GENOME_MAP[$SPECIES]}

    JXNHASH_INPUT="$PROJ/zenodo/datafiles/datafile_jxnHash_${GENOME}.csv"
    ERP_DATA_FILE="$PROJ/zenodo/datafiles/datafile_ERP_${GENOME}.csv"

    # Intermediates (documented names; deleted at end of iteration)
    UJC_PROPREADS="$ROZ/UJC_propReads_table_${SPECIES}.csv"
    ERP_PROPREADS="$ROZ/ERP_propReads_table_${SPECIES}.csv"
    UJC_FLAGGED="$ROZ/UJC_propReads_table_w_flags_${SPECIES}.csv"
    ERP_FLAGGED="$ROZ/ERP_propReads_table_w_flags_${SPECIES}.csv"
    UJC_GENE_COLS="$ROZ/gene_summary_from_UJC_${SPECIES}.csv"
    ERP_GENE_COLS="$ROZ/gene_summary_from_ERP_${SPECIES}.csv"
    GENE_SUMMARY="$ROZ/gene_summary_${SPECIES}.csv"
    GENE_SUMMARY_FLAGS=$ROZ/gene_summary_w_flags_${SPECIES}.csv

    # gene summary columns from AS analysis
    AS_GENE_SUM="$OUT/gene_summary_from_AS_analysis_${SPECIES}.csv"

    echo "Creating UJC propReads table"
    python "$SCRIPTS/create_UJC_propReads_table.py" \
        -i "$JXNHASH_INPUT" -o "$UJC_PROPREADS" -s "$SPECIES" -g "$GENOME" \
        -g1 "$GROUP1" -g2 "$GROUP2" -p "$COUNT_PREFIX" \
        --ujc-annotation-flag "$UJC_ANNOTATION_FLAG" --erp-annotation-flag "$ERP_ANNOTATION_FLAG"

    echo "Creating ERP propReads table"
    python "$SCRIPTS/create_ERP_propReads_table.py" \
        -i "$JXNHASH_INPUT" -o "$ERP_PROPREADS" -s "$SPECIES" \
        -g1 "$GROUP1" -g2 "$GROUP2" -p "$COUNT_PREFIX" --jxn-hash-col "${GENOME}_jxnHash" \
        --ujc-annotation-flag "$UJC_ANNOTATION_FLAG" --erp-annotation-flag "$ERP_ANNOTATION_FLAG"

    echo "Flagging UJC propReads table"
    python "$SCRIPTS/flag_annotation_status_and_bias_in_UJC_propReads_table.py" \
        -i "$UJC_PROPREADS" -o "$UJC_FLAGGED" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG/flag_UJC_propReads_${SPECIES}.log" \
        --ujc-annotation-flag "$UJC_ANNOTATION_FLAG" --erp-annotation-flag "$ERP_ANNOTATION_FLAG" \
        --group-class-col "$GROUP_CLASS_COL"

    echo "Flagging ERP propReads table"
    python "$SCRIPTS/flag_annotation_status_and_bias_in_ERP_propReads_table.py" \
        -i "$ERP_PROPREADS" -o "$ERP_FLAGGED" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG/flag_ERP_propReads_${SPECIES}.log" \
        --erp-data-file "$ERP_DATA_FILE" --group-class-col "$GROUP_CLASS_COL"

    echo "Aggregating UJC gene columns"
    python "$SCRIPTS/aggregate_UJC_propReads_to_gene_summary.py" \
        -i "$UJC_FLAGGED" -o "$UJC_GENE_COLS" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG/aggregate_UJC_gene_summary_${SPECIES}.log"

    echo "Aggregating ERP gene columns"
    python "$SCRIPTS/aggregate_ERP_propReads_to_gene_summary.py" \
        -i "$ERP_FLAGGED" -o "$ERP_GENE_COLS" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG/aggregate_ERP_gene_summary_${SPECIES}.log"

    echo "Detecting AS events (standalone table)"
    python "$SCRIPTS/detect_AS_events_from_UJC_propReads.py" \
        -i "$UJC_FLAGGED" -e "$ERP_FLAGGED" -o "$AS_GENE_SUM" -s "$SPECIES" \
        -g1 "$GROUP1" -g2 "$GROUP2" --group-class-col "$GROUP_CLASS_COL" --cpus "$CPUS"

    echo "Merging UJC and ERP_plus gene summary columns"
    python "$SCRIPTS/merge_gene_summary_columns.py" \
        --ujc-summary "$UJC_GENE_COLS" --erp-summary "$ERP_GENE_COLS" \
        -o "$GENE_SUMMARY" -s "$SPECIES" -l "$LOG/merge_gene_summary_columns_${SPECIES}.log"

    echo "Adding final flags and sex-bias categories"
    python "$SCRIPTS/add_flags_and_bias_categories_to_gene_summary.py" \
        -i "$GENE_SUMMARY" -o "$GENE_SUMMARY_FLAGS" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG/add_flags_gene_summary_${SPECIES}.log"

    echo "Completed $SPECIES"
done < <(tail -n +2 "$DESIGN")

echo ""
echo "Merging dmel6 ortholog geneID and gene symbol (all species)"
python "$SCRIPTS/merge_dmel6_geneSymbol_to_gene_summary.py" \
    --summary-dir "$OUT" \
    --in-prefix "gene_summary_w_flags" \
    --out-prefix "gene_summary_w_geneSymbol"

echo "Merging genesetID (all species)"
python "$SCRIPTS/merge_genesetID_to_gene_summaries.py" \
    --summary-dir "$OUT" \
    --in-prefix "gene_summary_geneSymbol" \
    --out-prefix "gene_summary"

echo ""
echo "  gene summary pipeline complete"
echo "Final outputs in: $OUT"
echo "  gene_summary_{species}.csv"
echo "  gene_summary_from_AS_analysis_{species}.csv"