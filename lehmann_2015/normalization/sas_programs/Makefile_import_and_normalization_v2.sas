libname CEGS '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname CEGLOCAL '!SASLOC1/cegs_sergey/sasdata';
libname DMEL '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Import Design File */
    * Since I am dealing waith a large number of datasets it will be easier if
    * I use a design file and iterate over them. The CEGS incomplete dataset
    * is a subset of the CEGS incomplete dataset.
    *
    * INPUT: !MCLAB/cegs_sergey/design_files/CEGS_70_lines_no_tech.txt
    *        !MCLAB/cegs_sergey/design_files/CEGS_incomplete_lines_no_tech.txt
    *
    * DATASET: CEGS.complete_design_by_rep
    *          CEGS.incomplete_design_by_rep
    *          CEGS.combined_design_by_rep
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/import_design_file.sas';

/* Flag Raleigh and Winters Lines */
    * It may be useful to have a simple flag_raleigh that I can merge on to
    * different datasets in the future.
    *
    * INPUT: CEGS.combined_design_by_rep
    *
    * DATASET CEGS.flag_raleigh
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_raleigh.sas';

/* Flag Dataset of origin */
    * There are several classifications of the CEGS data. (1) we have the
    * COMPLETE data that have at least 3 reps for each mating status. (2) There
    * is PARTIAL data that has 3 reps for at least one mating status. (3) There
    * is INCOMPLETE data that does not have 3 reps for a mating status.
    *
    * INPUT: CEGS.complete_design_by_rep
    *        CEGS.incomplete_design_by_rep
    *
    * DATASET CEGS.flag_dataset_of_origin
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_dataset_of_origin.sas';

/* Import Coverage Counts on Fusions */
    * Covearge counts on Fusions are imported iterating over the design
    * file and a final stacked is created. I have created coverage for the
    * complete dataset that we will be using for ASE.
    *
    * I have started using only reads that were unique, ie not duplicate sequence.
    *
    * INPUT: !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_non-redundant_fusions_on_fusions_nodup/&line._&mv.&rep..csv
    *        !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_non-redundant_fusions_on_fusions_incomplete_nodup/&line._&mv.&rep..csv
    *        CEGS.complete_design_by_rep
    *        CEGS.incomplete_design_by_rep
    *
    * DATASET: CEGS.ccfus_stack
    *          CEGLOCAL.ccfus_stack (LOCAL COPY FOR FASTER ACCESS)
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/import_coverage_counts_on_fusions.sas';

/* Import Coverage Counts on Junctions */
    * Covearge counts on Junctions are imported iterating over the design
    * file and a final stacked is created. I initially tried appending
    * junction counts to fusion counts, but this resulted in a very large
    * file (~30GB). I have decided to keep junction separate until the
    * normalization step.
    *
    * INPUT: !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_canonical_junctions_nodup/&line._&mv.&rep..csv
    *        !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_canonical_junctions_incomplete_nodup/&line._&mv.&rep..csv
    *        CEGS.complete_design_by_rep
    *        CEGS.incomplete_design_by_rep
    *
    * DATASET: CEGS.junc_cnts
    *          CEGLOCAL.junc_cnts (LOCAL COPY FOR FASTER ACCESS)
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/import_coverage_counts_on_junctions.sas';

/* Flag Sample On */
    * This script creates flags to show which line/fusion have expression. I
    * looked at several different measures, and we have decided to drop samples
    * that have <29300 exonic regions that have an APN > 0. These are low expressing samples.
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *
    * DATASET: CEGS.detected_exon_cnts
    *          CEGS.flag_sample_on
    *
    * FIGURE: !MLCAB/cegs_sergey/reports/line_normalization/exon_expression_counts.png
    *
    *                                                Cumulative    Cumulative
    * flag_sample_on    Frequency     Percent     Frequency      Percent
    * -------------------------------------------------------------------
    *              0          97       17.17            97        17.17
    *              1         468       82.83           565       100.00
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_sample_on_v2.sas';
    
/* Import contamination flags */
    * Using the ASE pipeline we have identified lines that appear to be
    * contaminated. I need to import these flags and create a permanent
    * dataset. 
    *
    * INPUT: !MCLAB/cegs_sergey/reports/check_ase_bias_cutoff100_check_list_final_flagged_v2.csv
    *        CEGS.combined_design_by_rep
    *
    * DATASET: CEGS.flag_contaminated
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_contaminated.sas';

/* Flag Fusions to drop from dataset prior to normalizing */
    * We have settled on UQ normalization, however prior to this we want to
    * drop any fusions that are not really expressed. 
    * 1) flag_fusion_on if APN > 0
    * 2) summaraize to fusion*line*mv if in >50% of reps then fusion was on.
    * 3) flag_drop_fusion if not on in 90% of line*mv
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGS.flag_sample_on
    *        CEGS.flag_contaminated
    *
    * DATASET: CEGS.flag_fusion_on
    *          CEGS.flag_drop_fusion
    *
    * MATED
    *                                              Cumulative    Cumulative
    * flag_drop_fusion    Frequency     Percent     Frequency      Percent
    * ---------------------------------------------------------------------
    *                0       31079       49.19         31079        49.19
    *                1       32102       50.81         63181       100.00
    *
    * VIRGIN
    *                                              Cumulative    Cumulative
    * flag_drop_fusion    Frequency     Percent     Frequency      Percent
    * ---------------------------------------------------------------------
    *                0       32722       51.79         32722        51.79
    *                1       30459       48.21         63181       100.00
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_fusion_on.sas';

/* Normalize Dataset on Fusions*/
    * Now that I have created large datasets with all of my sample information
    * I need to normalize (uq3 method). For now I am keeping the Mated and
    * Virgin Datasets together during the normalization process because I may
    * want to compare them in the future.
    *
    * NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 LINE 54!!!!!!!!!!
    * 
    *               mating_status=M         
    *
    *             The MEANS Procedure
    *
    *          Variable            Median
    *          --------------------------
    *          sum_mapped      4276504.01
    *          q3              33.8119518
    *          median          10.3346491
    *          --------------------------
    *
    *
    *               mating_status=V         
    *
    *          Variable            Median
    *          --------------------------
    *          sum_mapped      4263611.29
    *          q3              33.8645833
    *          median          10.8125000
    *          --------------------------
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGS.flag_sample_on
    *        CEGS.flag_drop_fusion
    *        CEGLOCAL.junc_cnts
    *
    * DATASET:  CEGS.norm_basic_stats_m
    *           CEGS.norm_basic_stats_v
    *           CEGLOCAL.norm_basic_stats_m
    *           CEGLOCAL.norm_basic_stats_v
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/normalize_log_uq3_v2.sas';

/* Check distributions */
    * Need to check the distribution of the data after normalization.
    *
    * INPUT: CEGLOCAL.norm_basic_stats_m 
    *        CEGLOCAL.norm_basic_stats_v
    *
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions_v2.R
    *
    * FILES: $MCLAB/cegs_sergey/reports/line_normalization/Mated_log_uq_apn_boxplot_v2.png
    *        $MCLAB/cegs_sergey/reports/line_normalization/Virgin_log_uq_apn_boxplot_v2.png
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/normalization_look_at_distributions.sas';
    
/* Merge on gene information */
    * Gene information is merged on by fusion_id. 
    *
    * INPUT:   CEGLOCAL.norm_basic_stats_m
    *          CEGLOCAL.norm_basic_stats_v
    *          DMEL.Fb551_si_fusions_unique_flagged
    *
    * DATASET: CEGLOCAL.ccfus_norm_stack_M 
    *          CEGLOCAL.ccfus_norm_stack_V
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/merge_gene_information_fusions.sas';
