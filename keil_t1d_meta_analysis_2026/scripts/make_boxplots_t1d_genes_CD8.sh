#!/bin/bash

PROJ=/nfshome/nkeil/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType
IND=${PROJ}/quantify_t1d_pacbio_transcripts
SCRIPTS=${PROJ}/scripts

GENE_LIST=${PROJ}/test_files_for_anna/t1d_gene_list_4_plots.csv
IN_CSV=${IND}/norm_data_frag_CD8_stack.csv
POSITIONS=${PROJ}/featureID_2_ef_start.csv
METRIC=log_uq_apn

OUTD=${PROJ}/boxplots_t1d_genes_CD8
    mkdir -p ${OUTD}

# Script paths
PER_FRAG_CC=$SCRIPTS/boxplot_expression_by_group_case_control_per_frag.py
PER_FRAG_MF=$SCRIPTS/boxplot_expression_by_group_case_control_M_F_per_frag.py
GENE_CC=$SCRIPTS/boxplot_expression_by_group_case_control.py
GENE_MF=$SCRIPTS/boxplot_expression_by_group_case_control_M_F.py


#Loop over genes
while IFS=, read -r SYMBOL ENSG _; do
  [[ -z "${SYMBOL:-}" || "$SYMBOL" =~ ^# ]] && continue
  [[ -z "${ENSG:-}" ]] && { echo "Skip line no ENSG for $SYMBOL"; continue; }

  GDIR=$OUTD/${ENSG}_${SYMBOL}
      mkdir -p "$GDIR"

  # 1) per-fragment case/control (PDF)
  python $PER_FRAG_CC \
      --in $IN_CSV \
      --gene $ENSG \
      --metric $METRIC \
      --out $GDIR/${ENSG}_${SYMBOL}_CD8_per_frag_case_control.pdf

  # 2) per-fragment M/F × case/control (PDF)
  python $PER_FRAG_MF \
      --in $IN_CSV \
      --gene $ENSG \
      --metric $METRIC \
      --out $GDIR/${ENSG}_${SYMBOL}_CD8_per_frag_M_F.pdf

  # 3) by-gene case/control (PNG)
  python $GENE_CC \
     --in $IN_CSV \
     --gene $ENSG \
     --metric $METRIC \
     --positions ${POSITIONS} \
     --out $GDIR/${ENSG}_${SYMBOL}_CD8_by_gene_case_control.png

  # 4) by-gene M/F × case/control (PNG)
  python $GENE_MF \
      --in $IN_CSV \
      --gene $ENSG \
      --metric $METRIC \
      --positions ${POSITIONS} \
      --out $GDIR/${ENSG}_${SYMBOL}_CD8_by_gene_M_F.png

done < $GENE_LIST

echo "Done"
