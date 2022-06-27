#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Make flags for chip rna supplemental file

## Get input directories and files
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms
IND=$PROJ/manuscript/supp_files

## Make output directory
OUTD=$PROJ/manuscript/results_for_ms_output
     mkdir -p ${OUTD}

## Tests for results in manuscript:
##     1) Faster X-hypothesis: more sex-biased expression on X vs. A
##        (Chi2 test using presence of at least 1 ttest significant
##        sex-biased fragment)

## NOTE: Was using features from $PROJ/results/data_files/mel_chip_rna_frag_flags_anno.csv
##       but these were missing many feautres includng all intergenic
python $PROJ/scripts/integrated_chip_rna/results_for_chip_rna_ms_04avn.py \
    -mg ${IND}/dmel_chip_rna_flags.csv \
    -mf $PROJ/ChIPseq/detection_above_background/features/mel_chip_5U_3U_TSS_frag_intr_inter_flag.csv \
    -sg ${IND}/dsim_chip_rna_flags.csv \
    -sf $PROJ/ChIPseq/detection_above_background/features/sim_chip_5U_3U_TSS_frag_intr_inter_flag.csv \
    -or ${IND}/dmel_dsim_ortholog_chip_rna_flags.csv \
    -d ${OUTD}
