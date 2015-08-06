libname CEGS '!MCLAB/cegs_sergey/sas_data';
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

/* Import Plate K Design File */
    * We have had an additional plate 'k' sequenced. I need to create a design
    * file for that plate and merge it the combined design file.
    *
    * INPUT: '!MCLAB/cegs_sergey/design_files/CEGS_platek_lines_no_tech.txt'
    *         CEGS.combined_design_by_rep
    *
    * DATASET: CEGS.platek_design
    *          CEGS.combined_w_platek_by_rep
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/import_platek_design_file.sas';

/* Flag Raleigh and Winters Lines */
    * It may be useful to have a simple flag_raleigh that I can merge on to
    * different datasets in the future.
    *
    * INPUT: CEGS.combined_w_platek_by_rep
    *
    * DATASET CEGS.flag_raleigh
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_raleigh_w_platek.sas';

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

/* Import Plate K */
    * Covearge counts on Fusions are imported iterating over the design
    * file and a final stacked is created. 
    *
    * I have started using only reads that were unique, ie not duplicate sequence.
    *
    * INPUT: !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_non-redundant_fusions_on_fusions_nodup/&line._&mv.&rep..csv
    *        !MCLAB/cegs_sergey/pipeline_output/OE_coverage_counts/coverage_count_fb551_non-redundant_fusions_on_fusions_incomplete_nodup/&line._&mv.&rep..csv
    *        CEGS.complete_design_by_rep
    *        CEGS.incomplete_design_by_rep
    *
    * DATASET: CEGS.platek_ccfus_stack
    *          CEGLOCAL.platek_ccfus_stack (LOCAL COPY FOR FASTER ACCESS)
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/import_coverage_counts_on_fusions_platek.sas';

/* Flag Sample On */
    * This script creates flags to show which line/fusion have expression. I
    * looked at several different measures, and we have decided to drop samples
    * that have <29300 exonic regions that have an APN > 0. These are low expressing samples.
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGLOCAL.platek_ccfus_stack 
    *
    * DATASET: CEGS.flag_sample_onk
    *
    *                                           Cumulative    Cumulative
    * flag_sample_on    Frequency     Percent     Frequency      Percent
    * -------------------------------------------------------------------
    *              0         104       16.46           104        16.46
    *              1         528       83.54           632       100.00
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_sample_on_w_platek.sas';
    
/* Import contamination flags */
    * Using the ASE pipeline we have identified lines that appear to be
    * contaminated. I need to import these flags and create a permanent
    * dataset. 
    *
    * INPUT: !MCLAB/cegs_sergey/reports/check_ase_bias_cutoff100_check_list_final_flagged_v2.csv
    *        CEGS.combined_design_by_rep
    *        CEGLOCAL.platek_ccfus_stack 
    *
    * DATASET: CEGS.flag_contaminatedk
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_contaminated_w_platek.sas';

/* Flag Fusions to drop from dataset prior to normalizing */
    * We have settled on UQ normalization, however prior to this we want to
    * drop any fusions that are not really expressed. 
    * 1) flag_fusion_on if APN > 0 after summing reps
    * 2) flag_drop_fusion_by_line if not on in 90% of line*mv
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGLOCAL.platek_ccfus_stack 
    *        CEGS.flag_sample_onk
    *        CEGS.flag_contaminatedk
    *
    * DATASET: CEGS.flag_drop_fusion_by_line
    *
    * MATED
    *                                              Cumulative    Cumulative
    * flag_drop_fusion    Frequency     Percent     Frequency      Percent
    * ---------------------------------------------------------------------
    *                0       37442       60.26         37442        60.26
    *                1       24689       39.74         62131       100.00
    *
    * VIRGIN
    *                                              Cumulative    Cumulative
    * flag_drop_fusion    Frequency     Percent     Frequency      Percent
    * ---------------------------------------------------------------------
    *                0       35676       57.49         35676        57.49
    *                1       26377       42.51         62052       100.00
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/flag_fusion_on_platek.sas';

/* Normalize Dataset on Fusions*/
    * Now that I have created large datasets with all of my sample information
    * I need to normalize (uq3 method). For now I am keeping the Mated and
    * Virgin Datasets together during the normalization process because I may
    * want to compare them in the future.
    *
    * NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 LINE 54!!!!!!!!!!
    * 
    *                mating_status=M           
    *
    *              The MEANS Procedure
    *
    *           Variable            Median
    *           --------------------------
    *           sum_mapped     11168302.12
    *           q3             156.9631579
    *           median          53.7710526
    *           --------------------------
    *
    *
    *                mating_status=V           
    *
    *           Variable            Median
    *           --------------------------
    *           sum_mapped      9193755.21
    *           q3             156.7578947
    *           median          52.7157895
    *           --------------------------
    *
    *
    * INPUT: CEGLOCAL.ccfus_stack
    *        CEGLOCAL.platek_ccfus_stack
    *        CEGS.flag_sample_onk
    *        CEGS.flag_drop_fusion_by_line
    *
    * DATASET:  CEGS.line_norm_basic_stats_m
    *           CEGS.line_norm_basic_stats_v
    *           CEGLOCAL.line_norm_basic_stats_m
    *           CEGLOCAL.line_norm_basic_stats_v
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/normalize_log_uq3_platek.sas';

/* Check distributions */
    * Need to check the distribution of the data after normalization.
    *
    * INPUT: CEGLOCAL.line_norm_basic_stats_m 
    *        CEGLOCAL.line_norm_basic_stats_v
    *
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions_v2.R
    *
    * FILES: $MCLAB/cegs_sergey/reports/line_normalization/Mated_log_uq_apn_boxplot_v2.png
    *        $MCLAB/cegs_sergey/reports/line_normalization/Virgin_log_uq_apn_boxplot_v2.png
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/line_normalization_look_at_distributions.sas';
