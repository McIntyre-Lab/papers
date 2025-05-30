libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Import parse motif results from MAST */
    * INPUT: '!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_a_results_up_and_down.csv'
    *        '!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_b_results_up_and_down.csv'
    *        '!MCLAB/arbeitman/arbeitman_Fru_network/motif_analysis/fru_c_results_up_and_down.csv'
    *
    * DATASET: FRU.motif_flags
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_import_motif_data_up_and_downstream.sas';

/* Calculate Motif Region Length */
    * Calculate the motif region length (2k+- gene start site);
    *
    * INPUT: DMEL530.fbgn2coord
    *
    * DATASET: FRU.motif_search_regions
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_create_gene_length_table.sas';

/* Look at Motif Region distributions with IND and REP */
    * There is some concern that we found induced genes having more motifs,
    * because induced genes may tend to be longer genes. It is expected to find
    * more motifs in longer genes by chance.
    *
    * Looks at the relationship between gene length and motif presence or
    * absence in a variety of ways.
    *
    * INPUT: FRU.motif_search_regions
    *        FRU.motif_flags_and_cnts
    *        FRU.Flag_ind_rep
    *
    * DATASET: FRU.motif_region_length_stats
    * 
    * FIGURES: '!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/length_figs/*'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_region_length_distributions.sas';

/* Check Gene Length Freqs with CMH */
    * Create contengency tables using CMH with Length as our strata value. 
    * 
    * INPUT: FRU.flag_ind_rep
    *        FRU.motif_flags_and_cnts
    *        FRU.motif_search_regions
    *
    * OUTPUTS a bunch of freqs
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_gene_length_check.sas';
    
/* Look at distribution of the number of motifs */
    * What is the distirubtion of Fru{A,B,C} motifs throughout the geneome. For
    * example, do genes have more FruA motifs then FruB.
    *
    * INPUT: FRU.motif_flags_and_cnts
    *
    * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_distribution.csv
    *
    * PLOTS: can be generated by->'!MCLAB/arbeitman/arbeitman_fru_network/r_programs/motif_distribtuion.R'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_motif_distribution.sas';


/* Male */
    /* Motif Analysis for Ind and Rep */
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed. Induced and
        * repressed were modeled separately.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_male
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_male.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_male.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_male.sas';

    /* Motif Analysis for Ind and Rep Version 2 */
        * Similar to above, except induced and repressed were modeled together.
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_male_v2
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_male_v2.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_male_v2.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_male_v2.sas';

    /* Multiple Motif Enrichment */
        * looked a little at motif number enrichment, but we decided that we
        * should not pursue this any further
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_multiple_male.sas';

    /* Look at distribution of the number of motifs */
        * Previously the global distribution of Fru{A,B,C} motifs was looked
        * at. Does the distribution look different for genes that are induced
        * or repressed in males?
        *
        * INPUT: FRU.motif_flags_and_cnts
        *        FRU.flag_ind_rep
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/male_&type._motif_distribution.csv
        *
        * PLOTS: can be generated by->'!MCLAB/arbeitman/arbeitman_fru_network/r_programs/motif_distribtuion.R'
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_male_motif_distribution.sas';


/* Female */

    /* Motif Analysis for Ind and Rep */
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed. Induced and
        * repressed were modeled separately.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_female
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_female.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_female.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_female.sas';

    /* Motif Analysis for Ind and Rep Version 2 */
        * Similar to above, except induced and repressed were modeled together.
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_female_v2
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_female_v2.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_female_v2.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_female_v2.sas';

    /* Multiple Motif Enrichment */
        * looked a little at motif number enrichment, but we decided that we
        * should not pursue this any futher
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_multiple_female.sas';

    /* Look at distribution of the number of motifs */
        * Previously the global distribution of Fru{A,B,C} motifs was looked
        * at. Does the distribution look different for genes that are induced
        * or repressed in females?
        *
        * INPUT: FRU.motif_flags_and_cnts
        *        FRU.flag_ind_rep
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/female_&type._motif_distribution.csv
        *
        * PLOTS: can be generated by->'!MCLAB/arbeitman/arbeitman_fru_network/r_programs/motif_distribtuion.R'
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_female_motif_distribution.sas';


/* Null */

    /* Motif Analysis for Ind and Rep */
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed. Induced and
        * repressed were modeled separately.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_null_male
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_null_male.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_null_male.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_null_male.sas';

    /* Motif Analysis for Ind and Rep Version 2 */
        * Similar to above, except induced and repressed were modeled together.
        * Fisher's exact test looking to see if there is enrichment of the
        * different binding motifs for genes induced or repressed.
        *
        * INPUT: FRU.flag_ind_rep
        *        FRU.motif_flags_and_cnts
        *        
        * DATASET: FRU.fru_motif_test_up_down_null_v2
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_&letter._motif_table_up_and_down_null_male_v2.csv
        *         !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_null_male_v2.csv
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_null_male_v2.sas';

    /* Multiple Motif Enrichment */
        * looked a little at motif number enrichment, but we decided that we
        * should not pursue this any futher
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_up_and_downstream_multiple_null_male.sas';

    /* Look at distribution of the number of motifs */
        * Previously the global distribution of Fru{A,B,C} motifs was looked
        * at. Does the distribution look different for genes that are induced
        * or repressed in null mutants?
        *
        * INPUT: FRU.motif_flags_and_cnts
        *        FRU.flag_ind_rep
        *
        * OUTPUT: !MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/null_male_&type._motif_distribution.csv
        *
        * PLOTS: can be generated by->'!MCLAB/arbeitman/arbeitman_fru_network/r_programs/motif_distribtuion.R'
        ;
        %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/motif_analysis_null_male_motif_distribution.sas';

