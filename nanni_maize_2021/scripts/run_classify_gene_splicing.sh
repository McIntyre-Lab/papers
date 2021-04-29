#!/bin/bash

## Must be run in env with matplotlib_venn installed
##     (used env called mylab on AVN grimshawi)

## Set directories
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio
DIST=$PROJ/sqanti_classification_category_subset/transcriptDistance_FSM_ISM_zmtr_NIC_NNC_monoexon_filter
EA=$PROJ/sqanti_classification_category_subset/EA_annotations
PLOT=$PROJ/sqanti_classification_category_subset/plot_consol_FSM_ISM_zmtr_NIC_NNC_monoexon_filter

## Classify types of splicing in genes
##     (alt. exon, alt. donor/acceptr, intron retention)
#python $PROJ/scripts/classify_gene_splicing_02avn.py \
#    -d ${DIST}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
#    -f ${EA}/zmtr_fsm_ism_nic_nnc_exon_fragment_annotations_flag_IR.csv \
#    -o ${DIST}/zmtr_fsm_ism_nic_nnc_gene_flag_splicing.csv


## Describe genes with 2 transcripts
#python $PROJ/scripts/describe_2_xcrpt_genes.py \
#    -f ${DIST}/zmtr_fsm_ism_nic_nnc_gene_flag_splicing.csv \
#    > ${DIST}/zmtr_fsm_ism_nic_nnc_describe_2_xcrpt_gene.txt

## Plot density for expression ratio of IR transcripts vs.
##     total expression within the gene in each genotype
python $PROJ/scripts/plot_IR_isoform_expression.py \
    -e $PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/flag_on_off_rsem_expression_isoforms_TPM_5_full_table.tsv \
    -d ${DIST}/zmtr_fsm_ism_nic_nnc_pairwise_transcript_distance_reduced_sharedJunction.csv \
    -f $PROJ/sqanti_classification_category_subset/EA_annotations/zmtr_fsm_ism_nic_nnc_exon_fragment_annotations_flag_IR.csv \
    -p ${DIST}/zmtr_fsm_ism_nic_nnc \
    > ${DIST}/zmtr_fsm_ism_nic_nnc_2xcrptGene_IR_counts.txt


## Get values for venn diagram
## Columns in output CSV file are:
## 1: 'gene_id'
## 2: 'flag_intron_retention'
## 3: 'flag_alt_exon'
## 4: 'flag_alt_donor_acceptor'
## 5: 'numXcrpt'

## Order of categories are Abc,aBc,ABc,abC,AbC,aBC,ABC
##     where A = alt exon, B = alt donor/acceptor, C = IR
#COUNTS=$(awk -F "," 'BEGIN{Abc=0;aBc=0;ABc=0;abC=0;AbC=0;aBC=0;ABC=0;} \
#                     NR!=1{if($3==1){ \
#                              if($4==1){ \
#                                 if($2==1){ABC=ABC+1} \
#                                 else{ABc=ABc+1}} \
#                              else{ \
#                                 if($2==1){AbC=AbC+AbC} \
#                                 else{Abc=Abc+1}}} \
#                           else{ \
#                              if($4==1){ \
#                                 if($2==1){aBC=aBC+1} \
#                                 else{aBc=aBc+1}} \
#                              else{ \
#                                 if($2==1){abC=abC+1} \
#                                 else{}}}} \
#                     END{print Abc","aBc","ABc","abC","AbC","aBC","ABC}' \
#             ${DIST}/zmtr_fsm_ism_nic_nnc_gene_flag_splicing.csv)
#echo ${COUNTS}

## Check counts add up porperly
#TOTAL=$(awk 'NR!=1' \
#    ${DIST}/zmtr_fsm_ism_nic_nnc_gene_flag_splicing.csv | wc -l)
#SUM=$(echo ${COUNTS} | awk -F "," '{print $1+$2+$3+$4+$5+$6+$7}')

#if [[ ${TOTAL} != ${SUM} ]]; then
#    echo "!!!ERROR: Counts do not match"
#fi

## Plot venn diagram
#python $PROJ/scripts/make_venn.py \
#    -n 3 \
#    -v "${COUNTS}" \
#    -l "Alt.Exons,Alt.Donor/Acceptor,IntronRetention" \
#    -o ${PLOT}/zmtr_fsm_ism_nic_nnc_gene_splicing_venn.png
