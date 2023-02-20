#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH
#export PATH=/Users/adalena/opt/anaconda3/envs/mylab/bin:$PATH

## Make flags for chip rna supplemental file

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
ORTHO=$PROJ/results/data_files/mel_sim_ortho_combo_flags_lit_gene_flag_GO_gene_flag.csv
MEL=$PROJ/RNAseq/model_output/mel_ttest_flags_with_anno_combo_flag_lit_gene_flag_GO_gene_flag.csv
SIM=$PROJ/RNAseq/model_output/sim_ttest_flags_with_anno_combo_flag.csv
CROSS=$PROJ/figures/mel_sim_ortho_map_scatter/ortho_map2bothCoord_all.csv
ANNOT=~/mclab/SHARE/McIntyre_Lab/useful_dmel_dsim_compare/dmel_orthologs_dsim_fb_2017_04.csv

## Make output directory
OUTD=$PROJ/supp_files

## Add gene-level files to supplementary gene-level result files
python $PROJ/scripts/integrated_chip_rna/add_chip_rna_paper_supp_flags.py \
    -m ${MEL} \
    -s ${SIM} \
    -or ${ORTHO} \
    -a ${ANNOT} \
    -x ${CROSS} \
    -d ${OUTD}

## Add feature files to supplementary feature-level result files
#python $PROJ/scripts/integrated_chip_rna/merge_features_4_supp.py \
#    -mr $PROJ/RNAseq/sample_summary/mel_5U_3U_TSS_frag_intr_inter_on_off_uq_ff.csv \
#    -sr $PROJ/RNAseq/sample_summary/sim_5U_3U_TSS_frag_intr_inter_on_off_uq_ff.csv \
#    -mc $PROJ/ChIPseq/detection_above_background/features/mel_chip_5U_3U_TSS_frag_intr_inter_flag.csv \
#    -sc $PROJ/ChIPseq/detection_above_background/features/sim_chip_5U_3U_TSS_frag_intr_inter_flag.csv \
#    -ma ~/mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/fb_features/mel_5U_3U_TSS_frag_intr_inter_uniq.csv \
#    -sa ~/mclab/SHARE/McIntyre_Lab/useful_dsim_data/flybase202/fb_features/sim_5U_3U_TSS_frag_intr_inter_uniq.csv \
#    -mp $PROJ/ChIPseq/detection_above_background/features/sample_level_DAI_files \
#    -sp $PROJ/ChIPseq/detection_above_background/features/sample_level_DAI_files \
#    -o ${OUTD}

