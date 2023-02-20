#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_ttest_flags_with_anno_combo_flag.csv
SIM=${INPUT}/sim_ttest_flags_with_anno_combo_flag.csv
MELFRAG=${INPUT}/mel_ttest_frags_with_anno.csv
SIMFRAG=${INPUT}/sim_ttest_frags_with_anno.csv


## Make directory for gene groups of interest
OUTPUT=${INPUT}/chip_rna_concordance_genes
    mkdir -p ${OUTPUT}
TSSOUT=${INPUT}/chip_rna_concordance_TSS
    mkdir -p ${TSSOUT}

for SPECIES in mel sim; do

    if [[ ${SPECIES} == "mel" ]]; then
        GENE=${MEL}
        FRAG=${MELFRAG}
    else
        GENE=${SIM}
        FRAG=${SIMFRAG}
    fi

    ## Concordance of gene sex bias using exonic regions
    ## Get frequencies of rna detected APN > 5 with gene_sex_bias
    ## gene_sex_bias is based off of APN > 0 detection with 2-fold change
    python ${SCRIPTS}/chip_rna_concordance_counts_05avn.py \
        -i ${GENE} \
        -p ${SPECIES} \
        -o ${INPUT}/${SPECIES}_chip_rna_concordance_counts.txt \
        -d ${OUTPUT}

    ## Concordance of TSS sex bias using unique TSS regions
#    python ${SCRIPTS}/chip_rna_concordance_TSS.py \
#        -i ${FRAG} \
#        -p ${SPECIES} \
#        -o ${INPUT}/${SPECIES}_chip_rna_concordance_TSS.txt \
#        -d ${TSSOUT}
done
