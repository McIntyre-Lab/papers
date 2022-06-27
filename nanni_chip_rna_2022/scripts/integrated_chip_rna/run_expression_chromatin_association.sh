#!/bin/sh

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/Dros_PB_ChIP
SCRIPTS=${PROJ}/scripts/integrated_chip_rna
INPUT=${PROJ}/RNAseq/model_output

## Get input files
MEL=${INPUT}/mel_ttest_flags_with_anno_combo_flag.csv
SIM=${INPUT}/sim_ttest_flags_with_anno_combo_flag.csv

for SPECIES in mel sim; do

    if [[ ${SPECIES} == "mel" ]]; then
        GENE=${MEL}
    else
        GENE=${SIM}
    fi

    ## Get associattions of sex-biased expression (split by X and A)
    ##     and the presence of male or female open chromatin
    ## Do Chi2 tests (test for association) and
    ##     Cramer's V (level of association)
    ## Check how many of the sex-biased expression is due to
    ##     association with the k27 of the opposite sex when there
    ##     is no k4 of the indicated sex
    python ${SCRIPTS}/expression_chromatin_association.py \
        -i ${GENE} \
        -o ${INPUT}/${SPECIES}_chip_rna_associations.txt
done

