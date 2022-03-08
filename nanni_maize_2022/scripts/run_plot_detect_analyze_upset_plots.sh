#!/bin/bash

## Must have matplotlib and UpSetPlot so using TranD conda env
export PATH=/TB14/TB14/conda_env/TranD/bin:$PATH

## Plot UpSet plots of detected and analyzable genes across genotypes
##     and output CSV file of flags and means
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL
SCRIPTS=$PROJ/2018/PacBio/scripts
IND=$PROJ/make_combination_flag_file

GTF=$PROJ/2018/PacBio/sqanti_classification_category_subset/ZMtr_from_pbID_fsm_ism_plus_nic_nnc_monoexon_filter_consol_junction.gtf

MEAN=$PROJ/2018/PacBio/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter/plot_maize_PB_groups_UpSet_v3/zmtr_fsm_ism_nic_nnc_consol_geneCollapse_from_meanDF.csv

CCS=$PROJ/2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_ccs_gene_htseq_count.tsv

SR=$PROJ/2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_shrtRd_gene_htseq_count.tsv

# Get flags
python ${SCRIPTS}/flag_detect_analyze_gene_file.py \
    -f ${IND}/combination_flag_file_shrtRd_ccs_hoopes.csv \
    -m ${MEAN} \
    -c ${CCS} \
    -s ${SR} \
    -g ${GTF} \
    -o ${IND}/maize_12604_flag_detect_analyze_mean.csv

# Plot UpSet plots
#python ${SCRIPTS}/plot_detect_analyze_upset_plots_03avn.py \
#    -i ${IND}/maize_12604_flag_detect_analyze_mean.csv \
#    -p ${IND}/maize_12604
