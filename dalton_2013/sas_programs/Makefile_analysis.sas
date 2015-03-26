libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Import coverage counts into SAS */
    * Import coverage counts and create a permanent sas dataset.
    *
    * INPUT: !MCLAB/arbeitman/arbeitman_fru_network/pipeline_output/coverage_on_fusions/all_coverage_counts.csv
    *
    * DATASET: FRU.all_coverage_counts      3617460 obs
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_import_coverage_counts_v2.sas';

/* Make design file */
    * Create a design file based upon an excel spreadsheet that was given to
    * use by Michelle. MCLAB/arbeitman/arbeitman_fru_network/documentation/Index Legend.xls
    *
    * DATASET: FRU.design_file
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_make_design_file_v2.sas';

/* Make stacked dataset */
    * Make a stacked dataset with log(RPKM). Any fusions that have a rpkm=0
    * will become undefined and will not be used in further modeling.
    *
    * INPUT: FRU.all_coverage_counts
    *        FRU.design_file
    *
    * DATASET: FRU.all_coverage_counts_with_key
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_make_stacked_with_key_v3.sas';

/* Flag No Treatment information */
    * Some of the data that I imported had not associated treatment
    * information. This was because we were only given treatment information
    * for the contrasts of interest. I filter this data of no interest here.  I
    * have added to this script so that it also creates a flag for treatments
    * with dsx in the treatment name. For now we are not wanting to include DSX
    * in the analysis so samples will be flagged and removed. I also noticed a
    * potential bug that would appear during the merge step of the foranalysis
    * dataset. I noticed that I was merging later on by only fusion_id. This is
    * not correct for flag_no_trt and flad_dsx because they are only uniq by
    * fusion_id*sample_id. So I made both flags to include fusion_id and
    * sample_id. 
    *
    *                                         Cumulative    Cumulative
    *   flag_no_trt    Frequency     Percent     Frequency      Percent
    *   ----------------------------------------------------------------
    *           0     2833677       78.33       2833677        78.33
    *           1      783783       21.67       3617460       100.00
    *
    *   783783/60291 = 13 samples without treatment information
    *
    * INPUT: FRU.all_coverage_counts_with_key
    *
    * DATASET: FRU.flag_no_trt
    *          FRU.flag_dsx
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_no_trt_v3.sas';

/* Drop Samples without Treatments and DSX from the dataset */
    * We are going to focus on FRU for this first paper, so we can go ahead and
    * drop samples without treatments and DSX.
    *
    * INPUT: FRU.all_coverage_counts_with_key
    *        FRU.flag_no_trt
    *        FRU.flag_dsx
    *
    * DATASET: WORK.short_list2
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_drop_no_trt_dsx.sas';

/* Flag No Variance */
    * We cannot analyze fusions that have no variance so I will go ahead and flag
    * these cases
    *
    *                                        Cumulative    Cumulative
    *   flag_no_var    Frequency     Percent     Frequency      Percent
    *   ----------------------------------------------------------------
    *          0       47401       78.62         47401        78.62
    *          1       12890       21.38         60291       100.00
    *
    *
    * INPUT: WORK.short_list2
    *
    * DATASET: FRU.flag_no_var
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_no_var_v3.sas';

/* Flag Zero Means */
    * This step is not needed I was just looking at the number of fusions that
    * had a mean of 0 RPKM
    * 
    *                                              Cumulative    Cumulative
    *   flag_zero_mean    Frequency     Percent     Frequency      Percent
    *   -------------------------------------------------------------------
    *                0       47401       78.62         47401        78.62
    *                1       12890       21.38         60291       100.00
    *
    * INPUT: WORK.short_list2
    *
    * DATASET: FRU.flag_zero_mean
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_zero_mean_v4.sas';

/* Make dataset for analysis */
    * Using my various flags I can not create a dataset for input into the
    * ANOVA model. I created a dataset with and without multigene fusions.
    *
    * INPUT: WORK.short_list2
    *        FRU.flag_no_var
    *        FRU.flag_zero_mean
    *        DMEL.fusion2gene2fbtr
    *
    * DATASET: FRU.foranalysis_no_multi
    *          FRU.foranalysis_with_multi
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_make_foranalysis_v3.sas';

/* ANOVA with Contrasts */
    * Run an ANOVA using log(RPKM).
    *
    * INPUT: FRU.foranalysis_no_multi
    *
    * DATASET: FRU.contrast_by_fusion
    *          FRU.lsmeans_by_fusion
    *          FRU.means_by_fusion
    *          FRU.model_fru_by_fusion
    *          FRU.resid_by_fusion
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_anova_v6.sas';

/* Flag Fusions Failing Normality */
    * Now we will look at our residuals and make sure they look ok before
    * proceeding.
    *
    *         flag_fail_                             Cumulative    Cumulative
    *          normality    Frequency     Percent     Frequency      Percent
    *   ---------------------------------------------------------------------
    *                  0       39419       83.16         39419        83.16
    *                  1        7982       16.84         47401       100.00
    *
    * INPUT: FRU.resid_by_fusion
    *
    * DATASET: FRU.flag_fail_normality_by_fusion
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_fail_normality_v2.sas';

/* Flag Missing */
    * Flag fusions that have missing data
    *
    * INPUT: FRU.resid_by_fusion
    *
    * DATASET: FRU.flag_missing
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_missing.sas';

/* ANOVA on gene level */
    * I tried running an ANOVA collapsed to the gene level, but there were too
    * many residuals that failed normality.
    *
    * INPUT: FRU.foranalysis_no_multi
    *
    * DATASET: FRU.contrast_by_symbol
    *          FRU.lsmeans_by_symbol
    *          FRU.means_by_symbol
    *          FRU.model_fru_by_symbol
    *          FRU.resid_by_symbol
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_anova_v7.sas';

/* Flag Genes Failing Normality */
    * Now we will look at our residuals and make sure they look ok before
    * proceeding. These look really bad so we won't proceed with this model.
    *
    *        flag_fail_                             Cumulative    Cumulative
    *          normality    Frequency     Percent     Frequency      Percent
    *   ---------------------------------------------------------------------
    *                  0        4211       40.84          4211        40.84
    *                  1        6101       59.16         10312       100.00
    *
    * INPUT: FRU.resid_by_gene
    *
    * DATASET: FRU.flag_fail_normality_by_gene
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_flag_fail_normality_v3.sas';

/* Calculate FDR for fusion Contrasts */
    * Next we will calculate the FDR on a fusion by fusion basis.
    *
    * INPUT: FRU.contrast_by_fusion
    *
    * DATASET: FRU.flag_fdr_contrasts_by_fusion
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_fdr_contrasts_v4.sas';

/* Calculate fusion level means */
    * Calculate the mean rpkm and log(rpkm) for all fusions.
    *
    * INPUT: FRU.all_coverage_counts_with_key
    *
    * DATASET: all_meansv2
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_means_v2';

/* Make Results Table */
    * Combine all results and flags to make final dataset.
    *
    * INPUT: DMEL.Fb530_si_fusions_unique_flagged
    *        FRU.model_fru_by_fusion
    *        FRU.all_meansv2
    *        FRU.flag_fdr_contrasts_by_fusion
    *        FRU.flag_fail_normality_by_fusion
    *        FRU.flag_no_var
    *        FRU.flag_zero_mean
    *        FRU.flag_missing
    *        DMEL.fusions2go
    *
    * DATASET: FRU.results_by_fusion
    *          FRU.results_plus_gov2
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_make_big_results_v6.sas';

/* Export Results Table */
    * Export results table to a CSV
    *
    * INPUT: FRU.results_plus_gov2
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/reports_external/gene_expression_results/results_plus_go_v2_20120721.csv'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_export_big_results_v5.sas';

/* Create Subset of Genes that failed normality */
    * There is a list of ~50 genes that are called out specifically in the
    * manuscript. I investigated this list of genes to see if any of them were
    * identified as being differentially expressed because of a fusion that was
    * also flagged as failing normality. I have identified 5 genes that meet
    * this criteria. Here I create a subset of these 5 genes to look at
    * further.
    * 
    * INPUT: FRU.foranalysis_no_multi
    *
    * DATASET: FRU.subset_fail_normality
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_create_subset_for_normality_checking.sas';

/* Look in more detail at specific fusions */
    * Now I will look in more detail at the residuals for specific fusions.
    * 
    * INPUT: FRU.subset_fail_normality
    *
    * DATASET:
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/analysis_investigate_fusions_fail_normality.sas';

