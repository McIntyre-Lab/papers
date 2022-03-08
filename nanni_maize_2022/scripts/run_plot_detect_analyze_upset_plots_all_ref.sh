#!/bin/bash

## Must have matplotlib and UpSetPlot so using TranD conda env
export PATH=/TB14/TB14/conda_env/TranD/bin:$PATH

## Plot UpSet plots of detected and analyzable genes across genotypes
##     and output CSV file of flags and means
PROJ=~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL
SCRIPTS=$PROJ/2018/PacBio/scripts
IND=$PROJ/make_combination_flag_file

#MEAN=$PROJ/2018/PacBio/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc/tappas_output_TMM_norm_1CPMfilter/plot_maize_PB_groups_UpSet_v3/zmtr_fsm_ism_nic_nnc_consol_geneCollapse_from_meanDF.csv

#CCS=$PROJ/2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_ccs_all_ref_gene_htseq_count.tsv

#SR=$PROJ/2018/PacBio/align_raw_reads_2_genomes/htseq_gene_counts/all_shrtRd_all_ref_gene_htseq_count.tsv

# Get flags
#python ${SCRIPTS}/flag_detect_analyze_gene_file.py \
#    -f ${IND}/combination_flag_file_shrtRd_ccs_hoopes.csv \
#    -m ${MEAN} \
#    -c ${CCS} \
#    -s ${SR} \
#    -o ${IND}/maize_b73_ref_flag_detect_analyze_mean.csv

# Plot UpSet plots
python ${SCRIPTS}/plot_detect_analyze_upset_plots_04avn.py \
    -i ${IND}/combination_flag_file_shrtRd_ccs_hoopes.csv \
    -d \
    -p ${IND}/maize_b73_ref
#python ${SCRIPTS}/plot_detect_analyze_upset_plots_03avn.py \
#    -i ${IND}/maize_b73_ref_flag_detect_analyze_mean.csv \
#    -p ${IND}/maize_b73_ref

## Get simple agreement between DE with rsem and DE with cvr
#python $PROJ/pacbio_paper/resubmission_2021/scripts/get_12604_rsem_vs_cvg_cnt_DE_counts.py

## Count number of detected/analyzable genes that are maize paralogs
## Columns:
##     1  geneID
##     5  flag_b73_reference_gene
##     6  flag_detect_ccs_b73_amb_gt0
##     7  flag_detect_ccs_b73_ele_gt0
##     8  flag_detect_ccs_b73_gt0
##     9  flag_detect_ccs_c123_amb_gt0
##    10  flag_detect_ccs_c123_ele_gt0
##    11  flag_detect_ccs_c123_gt0
##    12  flag_detect_ccs_hp301_amb_gt0
##    13  flag_detect_ccs_hp301_ele_gt0
##    14  flag_detect_ccs_hp301_gt0
##    15  flag_detect_ccs_mo17_amb_gt0
##    16  flag_detect_ccs_mo17_ele_gt0
##    17  flag_detect_ccs_mo17_gt0
##    18  flag_detect_ccs_nc338_amb_gt0
##    19  flag_detect_ccs_nc338_ele_gt0
##    20  flag_detect_ccs_nc338_gt0
##    21  flag_detect_shrtRd_b73_amb_gt0
##    22  flag_detect_shrtRd_b73_ele_gt0
##    23  flag_detect_shrtRd_b73_gt0
##    24  flag_detect_shrtRd_c123_amb_gt0
##    25  flag_detect_shrtRd_c123_ele_gt0
##    26  flag_detect_shrtRd_c123_gt0
##    27  flag_detect_shrtRd_hp301_amb_gt0
##    28  flag_detect_shrtRd_hp301_ele_gt0
##    29  flag_detect_shrtRd_hp301_gt0
##    30  flag_detect_shrtRd_mo17_amb_gt0
##    31  flag_detect_shrtRd_mo17_ele_gt0
##    32  flag_detect_shrtRd_mo17_gt0
##    33  flag_detect_shrtRd_nc338_amb_gt0
##    34  flag_detect_shrtRd_nc338_ele_gt0
##    35  flag_detect_shrtRd_nc338_gt0
##    38  flag_zea_mays_paralog_hoopes
##    52  flag_analyze_b73
##    53  flag_analyze_mo17
##    54  flag_analyze_c123
##    55  flag_analyze_hp301
##    56  flag_analyze_nc338
#awk -F "," '$5==1{}'  \
#    ${IND}/combination_flag_file_shrtRd_ccs_hoopes.csv
