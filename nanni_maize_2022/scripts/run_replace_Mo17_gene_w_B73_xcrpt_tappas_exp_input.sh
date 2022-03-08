#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Set directories
REF=~/mclab/SHARE/McIntyre_Lab/useful_maize_info
MO17=~/mclab/SHARE/McIntyre_Lab/useful_maize_mo17_data
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
IND=$PROJ/cvr_cnts_gene_mo17_cau

#GTF=${REF}/Zea_mays.B73_RefGen_v4.41.gtf.gz
#GTF=${REF}/Zea_mays.AGPv4.34.gtf.gz
G2T=${REF}/Zea_mays.AGPv4.34_gene_2_first_transcript.csv
SYN=${MO17}/synteny_Mo17Cau_B73v4/B73v4.36_Mo17CAU_synfind_synmap_SASoutput_avn.tsv

LOG=${IND}/replace_gene_w_transcript_tappas_v4.34.log
echo "#### Replace gene values with transcript values of B73 v4.34 ###" \
    > ${LOG}

for FILE in $(ls ${IND}/sbys_*_tpm.tsv); do
    PREFIX=$(dirname ${FILE})/$(basename ${FILE} .tsv)
    echo "
${FILE}" >> ${LOG}
    python $PROJ/scripts/replace_Mo17_gene_w_B73_xcrpt_tappas_exp_input.py \
        -e ${FILE} \
        -s ${SYN} \
        -g2t ${G2T} \
        -o ${PREFIX}_mod.tsv \
        >> ${LOG}
done
