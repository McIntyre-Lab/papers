#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL
IND=$PROJ/make_combination_flag_file/make_files_4_tappas_amm

#GTF=${REF}/Zea_mays.B73_RefGen_v4.41.gtf.gz
#GTF=${REF}/Zea_mays.AGPv4.34.gtf.gz
G2T=${REF}/Zea_mays.AGPv4.34_gene_2_first_transcript.csv

LOG=${IND}/replace_gene_w_transcript_tappas_v4.34.log
echo "#### Replace gene values with transcript values of B73 v4.34 ###" \
    > ${LOG}

for FILE in $(ls ${IND}/sbys_*_tpm.tsv); do
    PREFIX=$(dirname ${FILE})/$(basename ${FILE} .tsv)
    echo "
${FILE}" >> ${LOG}
    python ${IND}/replace_gene_w_transcript_tappas_exp_input.py \
        -e ${FILE} \
        -g2t ${G2T} \
        -o ${PREFIX}_mod.tsv \
        >> ${LOG}
done
