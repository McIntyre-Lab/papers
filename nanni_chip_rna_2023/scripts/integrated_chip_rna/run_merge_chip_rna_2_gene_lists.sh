#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Compare chip rna flags to gene lists

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
ORTHOTTEST=$PROJ/results/data_files/mel_sim_ortho_combo_flags.csv
MELTTEST=$PROJ/RNAseq/model_output/mel_ttest_flags_with_anno_combo_flag.csv
LISTDIR=~/mclab/SHARE/McIntyre_Lab/useful_dmel_data/gene_lists/dmel617_gene_lists
LITGENES=${LISTDIR}/dmel617_gene_list_flags.csv
GOGENES=${LISTDIR}/GO_associated_genes/dmel617_GO_gene_list_flags.csv
#SIM=$PROJ/CHIP_RNA_ms/results/data_files/sim_gene_flags_anno.csv

## Merge gene-level flags with gene lists from literature
python $PROJ/scripts/integrated_chip_rna/merge_chip_rna_2_gene_lists_03avn.py \
    -t ${MELTTEST} \
    -ot ${ORTHOTTEST} \
    -ml ${LITGENES} \
    -s lit_gene_flag

## Get output files from lit gene merge
TTESTlit=$(dirname ${MELTTEST})/$(basename ${MELTTEST} .csv)_lit_gene_flag.csv
ORTHOTTESTlit=$(dirname ${ORTHOTTEST})/$(basename ${ORTHOTTEST} .csv)_lit_gene_flag.csv

## Merge gene-level flags with gene lists from GO associations
python $PROJ/scripts/integrated_chip_rna/merge_chip_rna_2_gene_lists_03avn.py \
    -t ${TTESTlit} \
    -ot ${ORTHOTTESTlit} \
    -ml ${GOGENES} \
    -s GO_gene_flag

rm ${TTESTlit} ${ORTHOTTESTlit}
