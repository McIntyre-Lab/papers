#!/bin/bash
set -euo pipefail

# min10reads-only gene summary pipeline.
# Filters ro jxnHash with at least 10 reads summed across all samples
# then run full gene-summary pipeline.
# output:
#   gene_summary_min10reads_{species}.csv
#   gene_summary_min10reads_from_AS_analysis_{species}.csv

PROJ=/nfshome/mgaran/mclab/SHARE/McIntyre_Lab/sex_specific_splicing
SCRIPTS=$PROJ/scripts
OUT=$PROJ/zenodo/summary_files
LOG_DIR=$PROJ/Tables/log_files
DESIGN=$PROJ/design_files/fiveSpecies_M_F_designFile.csv
DATA_DIR=$PROJ/submission/supplementary_files/datafiles

UJC_ANNOTATION_FLAG=flag_jxnHash_in_fiveSpecies_full_anno
ERP_ANNOTATION_FLAG=flag_ERPanno
GROUP_NAME="Sex"
GROUP_CLASS_COL="sexClass"
CPUS=24
MIN_READS=10

PREFIX="gene_summary_min10reads"

declare -A GENOME_MAP=(
    [dmel]="dmel6"
    [dsim]="dsim2"
    [dser]="dser1"
    [dsan]="dsan1"
    [dyak]="dyak2"
)

echo "=========================================="
echo "min10reads gene summary pipeline"
echo "Group class column: $GROUP_CLASS_COL"
echo "min-reads threshold: $MIN_READS"
echo "Output prefix: $PREFIX"
echo "CPUs for AS detection: $CPUS"
echo "=========================================="

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

    JXNHASH_FULL="$PROJ/Tables/datafile_jxnHash_${GENOME}.csv"
    ERP_DATA_FILE="$DATA_DIR/datafile_erp_${GENOME}_w_annoFlag_flagNovel_02amm.csv"

    # Filtered jxnHash input (deleted at end of iteration)
    JXNHASH_INPUT="$OUT/datafile_jxnHash_${GENOME}_w_cols_min${MIN_READS}reads.csv"

    echo "Filtering jxnHash datafile to >= ${MIN_READS} reads"
    python "$SCRIPTS/filter_datafile_jxnHash_threshold_min10reads.py" \
        -i "$JXNHASH_FULL" \
        -o "$JXNHASH_INPUT" \
        -m "$MIN_READS" \
        -p "$COUNT_PREFIX"

    # Intermediates (documented names; deleted at end of iteration)
    UJC_PROPREADS="$OUT/UJC_propReads_table_${SPECIES}.csv"
    ERP_PROPREADS="$OUT/ERP_propReads_table_${SPECIES}.csv"
    UJC_FLAGGED="$OUT/UJC_propReads_table_w_flags_${SPECIES}.csv"
    ERP_FLAGGED="$OUT/ERP_propReads_table_w_flags_${SPECIES}.csv"
    UJC_GENE_COLS="$OUT/gene_summary_from_UJC_${SPECIES}.csv"
    ERP_GENE_COLS="$OUT/gene_summary_from_ERP_${SPECIES}.csv"
    GENE_SUMMARY="$OUT/gene_level_summary_${SPECIES}.csv"

    # Gene summary before symbol/genesetid merges (documented; consumed by the merge scripts)
    GENE_SUMMARY_FLAGGED="$OUT/gene_level_summary_w_flags_${SPECIES}.csv"

    # Kept final outputs (min10reads label carried via PREFIX)
    AS_GENE_COLS="$OUT/${PREFIX}_from_AS_analysis_${SPECIES}.csv"

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
        -l "$LOG_DIR/flag_UJC_propReads_${PREFIX}_${SPECIES}.log" \
        --ujc-annotation-flag "$UJC_ANNOTATION_FLAG" --erp-annotation-flag "$ERP_ANNOTATION_FLAG" \
        --group-class-col "$GROUP_CLASS_COL"

    echo "Flagging ERP propReads table"
    python "$SCRIPTS/flag_annotation_status_and_bias_in_ERP_propReads_table.py" \
        -i "$ERP_PROPREADS" -o "$ERP_FLAGGED" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG_DIR/flag_ERP_propReads_${PREFIX}_${SPECIES}.log" \
        --erp-data-file "$ERP_DATA_FILE" --group-class-col "$GROUP_CLASS_COL"

    echo "Aggregating UJC gene columns"
    python "$SCRIPTS/aggregate_UJC_propReads_to_gene_summary.py" \
        -i "$UJC_FLAGGED" -o "$UJC_GENE_COLS" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG_DIR/aggregate_UJC_gene_summary_${PREFIX}_${SPECIES}.log"

    echo "Aggregating ERP gene columns"
    python "$SCRIPTS/aggregate_ERP_propReads_to_gene_summary.py" \
        -i "$ERP_FLAGGED" -o "$ERP_GENE_COLS" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG_DIR/aggregate_ERP_gene_summary_${PREFIX}_${SPECIES}.log"

    echo "Detecting AS events (standalone table)"
    python "$SCRIPTS/detect_AS_events_from_UJC_propReads.py" \
        -i "$UJC_FLAGGED" -e "$ERP_FLAGGED" -o "$AS_GENE_COLS" -s "$SPECIES" \
        -g1 "$GROUP1" -g2 "$GROUP2" --group-class-col "$GROUP_CLASS_COL" --cpus "$CPUS"

    echo "Merging UJC and ERP_plus gene summary columns"
    python "$SCRIPTS/merge_gene_summary_columns.py" \
        --ujc-summary "$UJC_GENE_COLS" --erp-summary "$ERP_GENE_COLS" \
        -o "$GENE_SUMMARY" -s "$SPECIES" -l "$LOG_DIR/merge_gene_summary_columns_${PREFIX}_${SPECIES}.log"

    echo "Adding final flags and sex-bias categories"
    python "$SCRIPTS/add_flags_and_bias_categories_to_gene_summary_tables.py" \
        -i "$GENE_SUMMARY" -o "$GENE_SUMMARY_FLAGGED" -s "$SPECIES" -g1 "$GROUP1" -g2 "$GROUP2" \
        -l "$LOG_DIR/add_flags_gene_summary_${PREFIX}_${SPECIES}.log"

    echo "Cleaning up per-species intermediates for $SPECIES"
    rm -f \
        "$JXNHASH_INPUT" \
        "$UJC_PROPREADS" "$UJC_FLAGGED" "$ERP_PROPREADS" "$ERP_FLAGGED" \
        "$UJC_GENE_COLS" "$ERP_GENE_COLS" "$GENE_SUMMARY"

    echo "Completed $SPECIES"
done < <(tail -n +2 "$DESIGN")

echo ""
echo "Merging dmel6 ortholog geneID and gene symbol (all species)"
python "$SCRIPTS/merge_dmel6_geneSymbol_to_gene_summary.py" \
    --summary-dir "$OUT" \
    --in-prefix "gene_level_summary_w_flags" \
    --out-prefix "gene_level_summary_w_geneSymbol"

echo "Merging genesetID (all species)"
python "$SCRIPTS/merge_genesetID_to_gene_summaries.py" \
    --summary-dir "$OUT" \
    --in-prefix "gene_level_summary_w_geneSymbol" \
    --out-prefix "${PREFIX}"

echo "Cleaning up merge intermediates"
rm -f "$OUT/gene_level_summary_w_flags_"*.csv "$OUT/gene_level_summary_w_geneSymbol_"*.csv

echo ""
echo "min10reads gene summary pipeline complete"
echo "Final outputs in: $OUT"
echo "  gene_summary_min10reads_{species}.csv"
echo "  gene_summary_min10reads_from_AS_analysis_{species}.csv"