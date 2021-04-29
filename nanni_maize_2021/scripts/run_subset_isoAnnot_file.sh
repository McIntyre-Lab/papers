#!/bin/bash

export PATH=$HOME/conda/.envs/lab/bin:$PATH

## Get list of transcripts in consolidated maize PacBio transcriptome
##     (only transcripts detected in at least one sample with TPM>0 in 50% of replicates)
## Use previous tappas file and classification file to get the associated
##     transcripts for FSM/ISM isoforms in tappas file and replace with reference transcript_id
## Subset previous tappas file for list of detected transcripts in
##     consolidated maize PacBio transcriptome

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
SQANTI=$PROJ/sqanti_post_filter_b73
ANNOT=$PROJ/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc
OUTD=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
ROZ=${OUTD}/roz_subset_isoAnnot
    mkdir -p ${ROZ}

## Get list of transcripts in consolidated maize PacBio transcriptome
awk -F "\t" 'NR!=1{print $1}' ${OUTD}/sbys_B73_4_tappas_TPM.tsv \
    > ${ROZ}/transcript_list.txt

## Subset isoAnnot by list of transcripts
python $PROJ/scripts/subset_isoAnnot_file.py \
    -i ${ANNOT}/final_maize.annotation_output4TAPPAS_mod.gff \
    -c ${SQANTI}/SQANTI_classification.txt \
    -l ${ROZ}/transcript_list.txt \
    -o ${OUTD}/final_maize.annotation_output4TAPPAS_mod_consol.gff

rm -r ${ROZ}
