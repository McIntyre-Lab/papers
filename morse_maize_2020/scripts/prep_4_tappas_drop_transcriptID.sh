#!/bin/bash


MCLAB=/home/ammorse/mclab/SHARE/McIntyre_Lab
INPUT=$MCLAB/maize_ozone_FINAL/2018/PacBio/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc

for i in B73 C123 Mo17 Hp301 NC338
do
    sed -i 's/"//g' $INPUT/sbys_${i}_4_tappas.tsv
    sed -i 's/transcriptID/\t/g' $INPUT/sbys_${i}_4_tappas.tsv
done
