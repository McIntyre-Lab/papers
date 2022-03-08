#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Merge significant GO terms across enrichment tests

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
B73R=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter
B73CC=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/make_combination_flag_file/tappas_output
Mo17CC=$PROJ/tappas_output_mo17_cau_cvr_cnt
OUTD=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/pacbio_paper/genetics_2022/combined_GO_results
    mkdir -p ${OUTD}

## Merge sig GO terms
for TYPE in biological_process cellular_component molecular_function; do
    FILE1=${B73R}/GO_gse_tappas_DE_${TYPE}_sig_terms.csv
    NAME1=rsemB73map
    FILE2=${B73CC}/GO_gse_tappas_DE_${TYPE}_sig_terms.csv
    NAME2=ccB73map
    FILE3=${B73CC}/GO_gse_tappas_1to1BMo_DE_${TYPE}_sig_terms.csv
    NAME3=ccB73mapBMo1to1
    FILE4=${Mo17CC}/GO_gse_tappas_DE_${TYPE}_sig_terms.csv
    NAME4=ccMo17mapBMo1to1
    python $PROJ/scripts/merge_GO_sig_results.py \
        -f ${FILE1},${FILE2},${FILE3},${FILE4} \
        -n ${NAME1},${NAME2},${NAME3},${NAME4} \
        -p ${OUTD}/${TYPE} \
        -c ${NAME1},${NAME2} \
        -c ${NAME2},${NAME3} \
        -c ${NAME3},${NAME4} \
        -o ${OUTD}/comb_GO_gse_tappas_DE_${TYPE}_sig_terms.csv
done
