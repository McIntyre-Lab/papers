#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH
#export PATH=/Users/adalena/opt/anaconda3/envs/mylab/bin:$PATH

## Make flags for chip rna supplemental file

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
ORTHO=$PROJ/results/data_files/mel_sim_ortho_combo_flags_lit_gene_flag_GO_gene_flag.csv
MEL=$PROJ/RNAseq/model_output/mel_ttest_flags_with_anno_combo_flag_lit_gene_flag_GO_gene_flag.csv
SIM=$PROJ/RNAseq/model_output/sim_ttest_flags_with_anno_combo_flag.csv

## Make output directory
OUTD=$PROJ/manuscript/supp_files

## Add gene-level files to supplementary gene-level result files
python $PROJ/scripts/integrated_chip_rna/add_chip_rna_paper_supp_flags.py \
    -m ${MEL} \
    -s ${SIM} \
    -or ${ORTHO} \
    -d ${OUTD}
